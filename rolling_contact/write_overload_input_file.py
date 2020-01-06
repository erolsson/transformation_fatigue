from __future__ import print_function

import os

import numpy as np

from scipy.optimize import fmin
from scipy.special import ellipe, ellipk

from transformation_fatigue.input_file_reader.input_file_reader import InputFileReader
from transformation_fatigue.materials.materials import SS2506_no_trans
from transformation_fatigue.materials.materials import SS2506


def calculate_elliptic_eccentricity(R1, R2):
    if R1 < R2:
        R2, R1 = R1, R2

    def f(x):
        e = 1-x**2
        y = (R1/R2 - (ellipe(e)/x**2 - ellipk(e))/(ellipk(e) - ellipe(e)))**2
        return y
    return fmin(f, R2/R1)[0]


def calculate_elastic_contact_force(R1, R2, p0):
    x = calculate_elliptic_eccentricity(R1, R2)
    e = np.sqrt(1-x**2)
    F1 = (4/np.pi/e**2*x**1.5*((ellipe(e**2)/x**2-ellipk(e**2))*(ellipk(e**2) - ellipe(e**2)))**0.5)**(1./3)
    print(F1)
    E0 = 200E3/(1-0.3**2)
    R0 = np.sqrt(R1*R2)
    return F1**6*p0**3*np.pi**3*R0**2/6/E0**2


def create_roller_model(simulation_file_name, geometry_file_name, material, p0):
    reader = InputFileReader()
    reader.read_input_file(geometry_file_name)

    nodes_pos, elements_pos = reader.nodal_data, reader.elements
    nodes_neg = np.copy(nodes_pos)
    nodes_neg[:, 1] *= -1
    # Swapping the nodes in the z-direction to get correct ordering of the nodes
    elements_neg = {}
    for e_type, elements in elements_pos.items():
        element_data = np.copy(elements)
        temp_connectivity = element_data[:, 5:].copy()
        element_data[:, 5:] = element_data[:, 1:5]
        element_data[:, 1:5] = temp_connectivity
        elements_neg[e_type] = element_data

    reader.write_geom_include_file(simulation_directory + '/roller_part_x_pos.inc')
    surfaces = [('EXPOSED_SURFACE', 'EXPOSED_ELEMENTS'), ('X0_SURFACE', 'X0_ELEMENTS'), ('Z0_SURFACE', 'Z0_ELEMENTS'),
                ('INNER_SURFACE', 'INNER_ELEMENTS')]
    reader.write_sets_file(simulation_directory + '/roller_sets.inc', surfaces_from_element_sets=surfaces)
    reader.nodal_data = nodes_neg
    reader.elements = elements_neg
    reader.write_geom_include_file(simulation_directory + '/roller_part_x_neg.inc')

    file_lines = ['*Heading',
                  '\tModel of cylinder rolling on a rigid plane']
    for side in ['x_pos', 'x_neg']:
        lines = ['*Part, name=roller_' + side,
                 '\t*Include, Input=roller_part_' + side + '.inc',
                 '\t*Include, Input=roller_sets.inc',
                 '\t*Solid Section, elset=ALL_ELEMENTS, material=' + material.name,
                 '\t\t1.0',
                 '*End Part']
        file_lines.extend(lines)
    file_lines.append('*Part, name=rigid_plane')
    file_lines.append('*End Part')

    file_lines.extend(material.material_input_file_string())
    file_lines.append('*Initial Conditions, type=Solution, user')
    file_lines.append('*Initial Conditions, type=Stress, user')
    file_lines.append('*Initial conditions, type=temperature')
    file_lines.append('\tRoller_x_pos.ALL_NODES, 22.')
    file_lines.append('\tRoller_x_neg.ALL_NODES, 22.')
    file_lines.append('*Assembly, name=rolling_contact_model')
    overlap = 0.001
    d = 20.1 - overlap
    for side in ['x_pos', 'x_neg']:
        file_lines.append('\t*Instance, name=roller_' + side + ', part=roller_' + side)
        file_lines.append('\t\t0., 0., ' + str(d))
        file_lines.append('\t*End instance')
    file_lines.append('\t*Instance, name=rigid_plane, part=rigid_plane')
    file_lines.append('\t\t0., 0., 0.')
    file_lines.append('\t\t0., 0., 0., 1., 0., 0., 90')
    file_lines.append('\t\t*Node, nset=plane_ref_pt')
    file_lines.append('\t\t\t1, 0., 0., 0.')
    file_lines.append('\t\t*Surface, type=CYLINDER, name=rigid_plane')
    file_lines.append('\t\tSTART, -50.,  0.')
    file_lines.append('\t\tLINE,  50.,   0.')
    file_lines.append('\t\t*Rigid Body, ref node=plane_ref_pt, analytical surface=rigid_plane')
    file_lines.append('\t*End instance')

    file_lines.append('\t*Tie, name=tie_x0')
    file_lines.append('\t\troller_x_pos.x0_Surface, roller_x_neg.x0_Surface')
    file_lines.append('\t*Surface, name=coupling_surface, combine=union')
    file_lines.append('\t\troller_x_pos.z0_surface')
    file_lines.append('\t\troller_x_neg.z0_surface')
    file_lines.append('\t*Surface, name=contact_surface, combine=union')
    file_lines.append('\t\troller_x_pos.exposed_surface')
    file_lines.append('\t\troller_x_neg.exposed_surface')
    file_lines.append('\t*Node, nset=roller_ref_node')
    file_lines.append('\t\t900000, 0., 0., 20.')
    file_lines.append('\t*Coupling, Constraint name=roller_load_coupling, ref node=roller_ref_node, '
                      'surface=coupling_surface')
    file_lines.append('\t\t*Kinematic')
    file_lines.append('*End Assembly')
    file_lines.append('*Surface interaction, name=frictionless_contact')
    file_lines.append('*Contact pair, interaction=frictionless_contact, type=node to surface')
    file_lines.append('\tcontact_surface, rigid_plane.rigid_plane')
    file_lines.append('*Boundary')
    file_lines.append('\troller_ref_node, 1, 2')
    file_lines.append('\troller_ref_node, 4, 6')
    file_lines.append('\trigid_plane.plane_ref_pt, 1, 6')
    file_lines.append('\troller_x_pos.y0_nodes, 2, 2')
    file_lines.append('\troller_x_neg.y0_nodes, 2, 2')
    file_lines.append('*step, name=initiate_contact, nlgeom=Yes')
    file_lines.append('\t*Static')
    file_lines.append('\t\t1e-1, 1., 1e-12, 1.')
    force = calculate_elastic_contact_force(40.2/2, 46, p0)
    file_lines.append('\t*Cload')
    file_lines.append('\t\troller_ref_node, 3, ' + str(-force/2))
    file_lines.append('\t*Controls, reset')
    file_lines.append('\t*Controls, parameters=line search')
    file_lines.append('\t\t5, , , , ')
    file_lines.append('\t*Contact Controls, Stabilize')
    file_lines.append('\t*Contact Interference, shrink')
    file_lines.append('\t\tcontact_surface, rigid_plane.rigid_plane')
    file_lines.append('\t*Output, field')
    file_lines.append('\t\t*Element Output')
    file_lines.append('\t\t\tS, SDV')
    file_lines.append('\t\t*Node Output')
    file_lines.append('\t\t\tU, CF, RF, CM, RM')
    file_lines.append('\t\t*Contact  Output')
    file_lines.append('\t\t\tCSTRESS')
    file_lines.append('*End step')
    file_lines.append('*step, name=unloading, nlgeom=Yes')
    file_lines.append('\t*Static')
    file_lines.append('\t\t1e-2, 1., 1e-12, 1.')
    force = calculate_elastic_contact_force(40.2/2, 46, 0)
    file_lines.append('\t*Cload')
    file_lines.append('\t\troller_ref_node, 3, 0.')
    file_lines.append('\t*Output, field, time interval=0.1')
    file_lines.append('\t\t*Element Output')
    file_lines.append('\t\t\tS, SDV')
    file_lines.append('\t\t*Node Output')
    file_lines.append('\t\t\tU, CF, RF, CM, RM')
    file_lines.append('\t\t*Contact  Output')
    file_lines.append('\t\t\tCSTRESS')
    file_lines.append('*End step')

    with open(simulation_file_name, 'w') as input_file:
        for line in file_lines:
            input_file.write(line + '\n')
        input_file.write('**EOF')


if __name__ == '__main__':
    simulation_directory = os.path.expanduser('~/rolling_contact/mechanical_FEM/2_5GPa_overload_2GPa_nom/')
    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)

    model_file = os.path.expanduser('~/python_fatigue/rolling_contact/input_files/roller.inp')
    create_roller_model(simulation_directory + 'overload_no_trans.inp', model_file, SS2506_no_trans, 2500)
    create_roller_model(simulation_directory + 'overload.inp', model_file, SS2506, 2500)
