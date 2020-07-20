import os
from subprocess import Popen
import sys

try:
    import distro
except ImportError:
    import platform as distro

from transformation_fatigue.input_file_reader.input_file_reader import InputFileReader
from transformation_fatigue.materials.materials import SS2506

specimen = sys.argv[-2]
R = float(sys.argv[-1])

specimen_loads = {'smooth': {-1.: [737., 774., 820.], 0.: [425., 440.]},
                  'notched': {-1.: [427., 450.], 0.: [225., 240., 255.]}}


def write_mechanical_input_files(geom_include_file, directory, loads, no_steps=1, initial_inc=1e-2):
    input_file_reader = InputFileReader()
    input_file_reader.read_input_file(geom_include_file)
    input_file_reader.write_geom_include_file(directory + '/include_files/geom_pos.inc')
    input_file_reader.nodal_data[:, 2] *= -1
    load_nodes = input_file_reader.set_data['nset']['Specimen_load_nodes']
    support_nodes = input_file_reader.set_data['nset']['Specimen_support_nodes']
    x_sym_nodes = input_file_reader.set_data['nset']['Specimen_XSym_Nodes']
    y = -min([input_file_reader.nodal_data[n-1, 2] for n in x_sym_nodes])
    z = max([input_file_reader.nodal_data[n - 1, 3] for n in x_sym_nodes])

    load_pos = input_file_reader.nodal_data[load_nodes[0]-1, 1:4]
    support_pos = input_file_reader.nodal_data[support_nodes[0]-1, 1]
    wb = (2*y)**2*(2*z)/6

    for e_data in input_file_reader.elements.values():
        n = e_data.shape[1] - 1
        e_data[:, n//2 + 1: n + 1], e_data[:, 1: n//2 + 1] = e_data[:, 1: n//2 + 1], e_data[:, n//2 + 1: n + 1].copy()
    input_file_reader.write_geom_include_file(directory + '/include_files/geom_neg.inc')
    input_file_reader.write_sets_file(directory + '/include_files/set_data.inc',
                                      str_to_remove_from_setname='SPECIMEN_',
                                      surfaces_from_element_sets=[('YSYM_SURFACE', 'YSYM_ELEMENTS')])

    def write_inp_file(force):
        force_max = force*wb/(support_pos - load_pos[0])/2*(1 + (1 + R)/(1 - R))
        file_lines = ['*Heading',
                      '\tMechanical model for the fatigue specimen utmis_ ' + specimen]

        def write_part(part_sign):
            lines = ['*Part, name=' + 'specimen_part_' + part_sign,
                     '\t*Include, Input=include_files/geom_' + part_sign + '.inc',
                     '\t*Include, Input=include_files/set_data.inc',
                     '\t*Solid Section, elset=ALL_ELEMENTS, material=SS2506',
                     '\t\t1.0',
                     '*End Part']
            return lines

        file_lines += write_part('pos')
        file_lines += write_part('neg')

        file_lines.append('**')
        file_lines.append('*Amplitude, name=amp, time=total time')
        file_lines.append('\t0.0, 0.0')
        file_lines.append('\t1.0, 1.0')
        t = 1.
        for i in range(no_steps):
            file_lines.append('\t' + str(t + 1.0) + ', ' + str(float(R)))
            file_lines.append('\t' + str(t + 2.0) + ', 1.0')
            t += 2.
        file_lines.append('\t' + str(t + 1.) + ', 0.0')

        file_lines.append('*Assembly, name=pulsator_model')
        for sign in ['pos', 'neg']:
            file_lines.append('\t*Instance, name=specimen_part_' + sign + ' , part=specimen_part_' + sign)
            file_lines.append('\t*End Instance')

        file_lines.append('\t*Tie, name=y_plane')
        file_lines.append('\t\tspecimen_part_pos.ysym_surface, specimen_part_neg.ysym_surface')

        file_lines.append('\t*Node, nset=load_node')
        file_lines.append('\t\t999999, ' + str(load_pos[0]) + ', ' + str(-2*load_pos[1]) + ', 0.0')

        file_lines.append('\t*Surface, name=load_surface, Type=Node')
        file_lines.append('\t\tspecimen_part_pos.load_nodes')
        file_lines.append('\t*Coupling, Constraint name=load_node_coupling, '
                          'ref node=load_node, surface=load_surface')

        file_lines.append('\t\t*Kinematic')
        file_lines.append('\t\t2, 2')

        file_lines.append('*End Assembly')
        file_lines.extend(SS2506.material_input_file_string())
        for sign in ['pos', 'neg']:
            file_lines.append('*Boundary')
            file_lines.append('\tspecimen_part_' + sign + '.xsym_nodes,\tXSYMM')
            file_lines.append('\tspecimen_part_' + sign + '.zsym_nodes,\tZSYMM')

        file_lines.append('*Boundary')
        file_lines.append('\tspecimen_part_neg.support_nodes,\t2')

        file_lines.append('*Boundary')
        file_lines.append('\tload_node,\t1')
        file_lines.append('\tload_node,\t3, 6')

        file_lines.append('*Initial Conditions, type=Solution, user')
        file_lines.append('*Initial Conditions, type=Stress, user')
        file_lines.append('*Initial conditions, type=temperature')
        file_lines.append('\tspecimen_part_pos.ALL_NODES, 22')
        file_lines.append('\tspecimen_part_neg.ALL_NODES, 22')
        for step in range(no_steps):
            for direction in ['max_load', 'min_load']:
                step_name = 'step_' + str(step+1) + '_' + direction
                file_lines.append('*step, name=' + step_name + ', nlgeom=Yes, inc=100000')
                file_lines.append('\t*Static')
                file_lines.append('\t\t' + str(initial_inc) + ', 1., 1e-12, 1.')
                file_lines.append('\t*CLoad, Amplitude=amp')
                file_lines.append('\t\tload_node, 2, ' + str(-force_max))
                file_lines.append('\t*Output, field')
                file_lines.append('\t\t*Element Output')
                file_lines.append('\t\t\tS, SDV, PEEQ')
                file_lines.append('\t\t*Node Output')
                file_lines.append('\t\t\tU')
                file_lines.append('*End step')

        file_lines.append('*step, name=relax, nlgeom=Yes, inc=100000')
        file_lines.append('\t*Static')
        file_lines.append('\t\t' + str(initial_inc) + ', 1., 1e-12, 1.')
        file_lines.append('\t*CLoad, Amplitude=amp')
        file_lines.append('\t\tload_node, 2, 0.')
        file_lines.append('\t*Output, field')
        file_lines.append('\t\t*Element Output')
        file_lines.append('\t\t\tS, SDV, PEEQ')
        file_lines.append('\t\t*Node Output')
        file_lines.append('\t\t\tU')
        file_lines.append('*End step')

        job_name = 'utmis_' + specimen + '_' + str(load).replace('.', '_') + '_R=' + str(int(R))
        with open(directory + '/' + job_name + '.inp', 'w') as input_file:
            for line in file_lines:
                input_file.write(line + '\n')
        return job_name
    job_names = []
    for load in loads:
        job_names.append(write_inp_file(load))
    return job_names


def write_run_file(job_names, heat_treatment_file, directory, cpus=12):
    file_lines = ['#!/bin/bash']
    if distro.linux_distribution()[0] == 'Ubuntu':
        abq = '\"singularity exec --nv ' + os.path.expanduser('~/imgs/sing/abaqus-2018-centos-7.img') + \
                   ' vglrun /opt/abaqus/2018/Commands/abq2018\"'
    else:
        abq = '/scratch/users/erik/SIMULIA/CAE/2018/linux_a64/code/bin/ABQLauncher'
        file_lines.extend(['#PBS -V',
                           '#PBS -z',
                           '#PBS -l select=1:ncpus=' + str(cpus),
                           'cd $PBS_O_WORKDIR'])

    file_lines.append('abq=' + abq)
    file_lines.append('subroutine=' + os.path.expanduser('~/python_projects/transformation_fatigue/'
                                                         'transformation_subroutine/transformation_subroutine.o'))

    for job_name in job_names:
        file_lines.append('cp ' + heat_treatment_file + ' ' + job_name + '.htd')
        file_lines.append('${abq} j=' + job_name + ' cpus=' + str(cpus) + ' user=${subroutine} interactive')
    with open(directory + '/run_utmis_' + specimen + '_R=' + str(int(R)) + '.sh', 'w') as shell_file:
        for line in file_lines:
            shell_file.write(line + '\n')


if __name__ == '__main__':
    heat_treatment_data_file = os.path.expanduser('~/utmis_specimens/smooth/Toolbox_Cooling_utmis_smooth.htd')
    simulation_directory = os.path.expanduser('~/utmis_specimens/smooth/mechanical_analysis/')
    geom_filename = os.path.expanduser('~/python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                       + specimen + '/utmis_' + specimen + '.inc')

    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)
    if not os.path.isdir(simulation_directory + '/include_files'):
        os.makedirs(simulation_directory + '/include_files')

    specimen_name = 'utmis_' + specimen

    jobs = write_mechanical_input_files(geom_filename, simulation_directory, specimen_loads[specimen][R], no_steps=2)
    write_run_file(job_names=jobs, directory=simulation_directory, heat_treatment_file=heat_treatment_data_file)
    # current_directory = os.getcwd()
    # os.chdir(simulation_directory)
    # Popen('qsub run_utmis_' + specimen + '_R=' + str(int(R)) + '.sh', shell=True)
    # os.chdir(current_directory)
