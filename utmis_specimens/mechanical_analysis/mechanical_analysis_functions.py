from collections import namedtuple
import os
import pathlib

try:
    import distro
except ImportError:
    import platform as distro

import numpy as np

from input_file_reader.input_file_reader import InputFileReader
from input_file_reader.input_file_functions import mirror_model


class Step:
    def __init__(self, name, load, max_inc=1., output_frequency=1):
        self.name = name
        self.load = load
        self.max_inc = max_inc
        self.output_frequency = output_frequency


Simulation = namedtuple('Simulation', ['name', 'steps', 'mode'])


def write_mechanical_input_files(specimen, geom_include_file, directory, simulations, material, initial_inc=1e-2):
    input_file_reader = InputFileReader()
    input_file_reader.read_input_file(geom_include_file)

    x = input_file_reader.nodal_data[:, 1]

    r = np.sqrt(np.sum(input_file_reader.nodal_data[:, 1:4]**2, 1))
    clamped_nodes = input_file_reader.nodal_data[x > 15, 0]
    center_node = input_file_reader.nodal_data[r < 1e-4, 0]

    for set_type in ['nset', 'elset']:
        for name in list(input_file_reader.set_data[set_type].keys()):
            new_name = name.lower().replace('specimen_', '')
            input_file_reader.set_data[set_type][new_name] = input_file_reader.set_data[set_type].pop(name)
    exposed_nodes = input_file_reader.set_data['nset']['exposed_nodes']
    y_sym_nodes = input_file_reader.set_data['nset']['ysym_nodes']
    clamped_nodes = list((set(clamped_nodes) & set(exposed_nodes)) - set(y_sym_nodes))
    input_file_reader.create_node_set('clamped_nodes', clamped_nodes)
    input_file_reader.create_node_set('center_node', [center_node])
    if not (directory / 'include_files').is_dir():
        (directory / 'include_files').mkdir(parents=True)
    input_file_reader.write_geom_include_file(directory / 'include_files/geom_pos.inc')
    input_file_reader.write_sets_file(directory / 'include_files/set_data.inc',
                                      str_to_remove_from_setname='SPECIMEN_',
                                      surfaces_from_element_sets=['ysym'])
    mirror_model(input_file_reader, 'y')
    input_file_reader.write_geom_include_file(directory / 'include_files/geom_neg.inc')

    def write_inp_file(load_point, simulation_data):
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

        file_lines.append('*Assembly, name=pulsator_model')
        for sign in ['pos', 'neg']:
            file_lines.append('\t*Instance, name=specimen_part_' + sign + ' , part=specimen_part_' + sign)
            file_lines.append('\t*End Instance')

        file_lines.append('\t*Tie, name=y_plane')
        file_lines.append('\t\tspecimen_part_pos.ysym_surface, specimen_part_neg.ysym_surface')

        file_lines.append('\t*Node, nset=load_point')
        file_lines.append('\t\t999999, ' + str(load_point[0]) + ', ' + str(load_point[1]) + ', ' + str(load_point[2]))

        file_lines.append('\t*Node, nset=symmetry_point')
        file_lines.append('\t\t999998, ' + str(load_point[0]) + ', ' + str(load_point[1]) + ', ' + str(load_point[2]))

        file_lines.append('\t*Surface, name=clamped_surface_pos, Type=Node')
        file_lines.append('\t\tspecimen_part_pos.clamped_nodes')
        file_lines.append('\t*Surface, name=clamped_surface_neg, Type=Node')
        file_lines.append('\t\tspecimen_part_neg.clamped_nodes')
        file_lines.append('\t*Surface, name=clamped_surface, Combine=Union')
        file_lines.append('\t\tclamped_surface_pos, clamped_surface_neg')

        file_lines.append('\t*Surface, name=xsym_surface_pos, Type=Node')
        file_lines.append('\t\tspecimen_part_pos.xsym_nodes')
        file_lines.append('\t*Surface, name=xsym_surface_neg, Type=Node')
        file_lines.append('\t\tspecimen_part_neg.xsym_nodes')
        file_lines.append('\t*Surface, name=xsym_surface, Combine=Union')
        file_lines.append('\t\txsym_surface_pos, xsym_surface_neg')

        file_lines.append('\t*Coupling, Constraint name=load_node_coupling, '
                          'ref node=load_point, surface=clamped_surface')
        file_lines.append('\t\t*Kinematic')
        file_lines.append('\t\t1, 6')

        file_lines.append('\t*Coupling, Constraint name=load_node_coupling, '
                          'ref node=symmetry_point, surface=xsym_surface')
        file_lines.append('\t\t*Kinematic')
        file_lines.append('\t\t1, 6')

        file_lines.append('*End Assembly')
        file_lines.extend(material.material_input_file_string())
        for sign in ['pos', 'neg']:
            file_lines.append('*Boundary')
            file_lines.append('\tspecimen_part_' + sign + '.zsym_nodes,\tZSYMM')

        file_lines.append('*Boundary')
        file_lines.append('\tload_point, 1, 5')
        file_lines.append('\tsymmetry_point, 1, 3, 4, 5')
        file_lines.append('*Initial Conditions, type=Solution, user')
        file_lines.append('*Initial Conditions, type=Stress, user')
        file_lines.append('*Initial conditions, type=temperature')
        file_lines.append('\tspecimen_part_pos.ALL_NODES, 22')
        file_lines.append('\tspecimen_part_neg.ALL_NODES, 22')
        for step in simulation_data.steps:
            file_lines.append('*step, name=' + step.name + ', nlgeom=Yes, inc=100000')
            file_lines.append('\t*Static')
            file_lines.append('\t\t' + str(initial_inc) + ', 1., 1e-12, ' + str(step.max_inc))
            if simulation.mode == 'displacement':
                file_lines.append('\t*Boundary')
                file_lines.append('\t\tload_point, 6, 6, ' + str(step.load))
            elif simulation.mode == 'force':
                file_lines.append('\t*Cload')
                file_lines.append('\t\tload_point, 6, ' + str(step.load))
            file_lines.append('\t*Output, field, frequency=' + str(step.output_frequency))
            file_lines.append('\t\t*Element Output')
            file_lines.append('\t\t\tS, SDV')
            file_lines.append('\t\t*Node Output')
            file_lines.append('\t\t\tU, RF, CF')
            file_lines.append('*End step')

        job_name = 'utmis_' + specimen + '_' + simulation_data.name
        with open(directory / (job_name + '.inp'), 'w') as input_file:
            for line in file_lines:
                input_file.write(line + '\n')
        return job_name
    job_names = []
    for simulation in simulations:
        job_names.append(write_inp_file([0, 0, 0], simulation))
    return job_names


def write_run_file(job_names, heat_treatment_file, run_file_name, cpus=12):
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
    file_lines.append('subroutine=' + str(pathlib.Path.home() / 'python_projects/transformation_fatigue'
                      / 'transformation_subroutine/transformation_subroutine.o'))

    for job_name in job_names:
        file_lines.append('cp ' + str(heat_treatment_file) + ' ' + job_name + '.htd')
        file_lines.append('${abq} j=' + job_name + ' cpus=' + str(cpus) + ' user=${subroutine} interactive')
    with open(run_file_name, 'w') as shell_file:
        for line in file_lines:
            shell_file.write(line + '\n')
