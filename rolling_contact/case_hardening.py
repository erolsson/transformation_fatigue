from collections import namedtuple
import os
import shutil

from case_hardening_toolbox import CaseHardeningToolbox
from case_hardening_toolbox import write_geometry_files_for_dante

from diffusivity import write_diffusion_file

from case_hardening_materials.gear_materials import SS2506

"""
    This file can be used as a template for setting up heat treatment simulations. The file creates one directory for 
    each simulation defined in the list simulations. Common files are placed in the directory include_file_directory 

    For setting up a simulation an input file is needed given as the first argument to write_geometry_files_for_dante

    The input file must contain:
        All needed nodes
        All needed elements
        A node set EXPOSED_NODES that are surface nodes in the carburization process
        An element set EXPOSED_SURFACE that are the elements corresponding to EXPOSED_NODES. These elements are used 
        for defining a surface corresponding to EXPOSED_NODES

        Node sets needed for boundary conditions. The boundary conditions should be specified in an own file and is 
        copied to the include_file_directory on lines 54-55
"""

Simulation = namedtuple('Simulation', ['simulation_directory', 'times', 'temperatures', 'carbon', 'tempering'])

# This is the main directory where all simulation folders will be placed
base_directory = os.path.expanduser('~/rolling_contact/')

name = 'roller'

simulations = [Simulation(simulation_directory='1',
                          times=[300., 90.], temperatures=[930, 850], carbon=[1.1, 0.75], tempering=(170, 120))]

# In this directory all common files for all heat treatment simulations will be placed
include_file_directory = base_directory + '/include_files/'

# This file contains all nodes, elements and sets to be further processed by the script
geometry_file_name = os.path.expanduser('~/python_fatigue/rolling_contact/input_files/roller.inp')

# All included files will be named after this like include_file_name_geo.inc or include_file_name_sets.inc
include_file_name = name

# Path to the boundary conditions file, copied to the include_file_directory
bc_file = os.path.expanduser('~/python_fatigue/rolling_contact/input_files/roller_bc.inc')

# Path to interaction property file
interaction_property_file = 'data_files/interaction_properties.inc'

# Name of the diffusion file, will be created in include_file_directory
diffusion_file_name = 'diffusivity.inc'

# Material, needs a python dict with material composition to generate the diffusion file
material = SS2506.composition

if not os.path.isdir(base_directory):
    os.makedirs(base_directory)

if not os.path.isdir(include_file_directory):
    os.makedirs(include_file_directory)

shutil.copyfile(bc_file, include_file_directory + '/' + include_file_name + '_BC.inc')

# writes the necessary geometry files and set files to include_file_directory
write_geometry_files_for_dante(geometry_data_file=geometry_file_name,
                               directory_to_write=include_file_directory,
                               dante_include_file_name=include_file_name,
                               boundary_condition_file=bc_file,
                               str_to_remove_from_set_names='Specimen_')

# Copying the interaction property file to the include file directory
shutil.copyfile(interaction_property_file, include_file_directory + '/interaction_properties.inc')

for simulation in simulations:
    inc_file_directory = os.path.relpath(include_file_directory, base_directory + '/' + simulation.simulation_directory)
    toolbox_writer = CaseHardeningToolbox(name=name,
                                          include_file_name=include_file_name,
                                          include_file_directory=inc_file_directory, env='KTH')
    toolbox_writer.cpus = 8
    toolbox_writer.diffusion_file = 'diffusivity.inc'
    toolbox_writer.interaction_property_file = 'interaction_properties.inc'
    toolbox_writer.heating_data.carbon = 0.5
    toolbox_writer.heating_data.time = 90.
    toolbox_writer.heating_data.temperature = 930.

    toolbox_writer.quenching_data.time = 3600.
    toolbox_writer.quenching_data.temperature = 120.

    toolbox_writer.cooldown_data.temperature = 80
    toolbox_writer.cooldown_data.time = 3600

    toolbox_writer.material = 'U925062'

    toolbox_writer.tempering_data.temperature = simulation.tempering[0]
    toolbox_writer.tempering_data.time = simulation.tempering[1]

    toolbox_writer.max_temp_inc = 5.
    toolbox_writer.max_vf_inc = 0.05
    toolbox_writer.write_dctrl_file = True

    toolbox_writer.add_carburization_steps(times=simulation.times, temperatures=simulation.temperatures,
                                           carbon_levels=simulation.carbon)
    directory_name = base_directory + '/' + simulation.simulation_directory
    current_directory = os.getcwd()
    if not os.path.isdir(directory_name):
        os.makedirs(directory_name)
    os.chdir(directory_name)
    toolbox_writer.write_files()
    os.chdir(current_directory)

write_diffusion_file(filename=include_file_directory + diffusion_file_name,
                     material=material)
