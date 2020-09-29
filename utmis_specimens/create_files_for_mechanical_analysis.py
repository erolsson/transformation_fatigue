from __future__ import division, print_function
import os
from collections import namedtuple
from subprocess import Popen
import sys

try:
    import distro
except ImportError:
    import platform as distro

import numpy as np

from transformation_fatigue.input_file_reader.input_file_reader import InputFileReader
from transformation_fatigue.materials.materials import SS2506
# from transformation_fatigue.materials.materials import SS2506Elastic as SS2506

# specimen = sys.argv[-2]
# R = float(sys.argv[-1])
specimen = 'smooth'
R = -1

specimen_loads = {'smooth': {-1.: [737., 774., 820.], 0.: [425., 440.]},
                  'notched': {-1.: [427., 450.], 0.: [225., 240., 255.]}}







def main():
    heat_treatment_data_file = os.path.expanduser('~/utmis_specimens/' + specimen + '/Toolbox_Cooling_utmis_' + specimen
                                                  + '.htd')
    simulation_directory = os.path.expanduser('~/utmis_specimens/' + specimen + '/mechanical_analysis_no_trans/')
    geom_filename = os.path.expanduser('~/python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                       + specimen + '/utmis_' + specimen + '.inc')

    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)
    if not os.path.isdir(simulation_directory + '/include_files'):
        os.makedirs(simulation_directory + '/include_files')

    specimen_name = 'utmis_' + specimen
    smax = 1000.
    Wb = 2.*5.**2/6.
    simulations = [Simulation(name='loading', steps=[Step('loading', smax*Wb)], mode='force')]
    jobs = write_mechanical_input_files(geom_filename, simulation_directory, simulations)
    write_run_file(job_names=jobs, directory=simulation_directory, heat_treatment_file=heat_treatment_data_file)


if __name__ == '__main__':
    main()
