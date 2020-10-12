import os
import sys

import numpy as np

from mechanical_analysis_functions import write_mechanical_input_files, Simulation, Step, write_run_file
from transformation_fatigue.materials.materials import SS2506


def main():
    wb = 2*25/6
    no_steps = 1
    specimen = sys.argv[-2]
    R = float(sys.argv[-1])
    specimen_loads = {'smooth': {-1.: [737., 774., 820.], 0.: [425., 440.]},
                      'notched': {-1.: [427., 450.], 0.: [225., 240., 255.]}}

    simulations = []
    for load_amplitude in specimen_loads[specimen][R]:
        mean_load = (1 + R)/(1 - R)*load_amplitude
        max_load = mean_load + load_amplitude
        min_load = mean_load - load_amplitude

        steps = []
        for step in range(1, no_steps + 1):
            steps.append(Step(str(step) + "_max_load", max_load*wb, output_frequency=1.))
            steps.append(Step(str(step) + "_min_load", min_load*wb, output_frequency=1.))
        steps.append(Step("relax", 0., output_frequency=1.))
        simulations.append(Simulation("snom=" + str(int(load_amplitude)) + "_R=" + str(int(R)), steps, 'force'))
    geom_filename = os.path.expanduser('~/python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                       + specimen + '/utmis_' + specimen + '.inc')
    simulation_directory = os.path.expanduser('~/utmis_specimens/' + specimen
                                              + '/mechanical_analysis/force_control_2/')
    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)
    job_names = write_mechanical_input_files(specimen, geom_filename, simulation_directory, simulations, SS2506)
    heat_treatment_data_file = os.path.expanduser('~/utmis_specimens/' + specimen + '/Toolbox_Cooling_utmis_'
                                                  + specimen + '.htd')
    write_run_file(job_names, heat_treatment_data_file, simulation_directory + '/run_R=' + str(int(R)) + '.sh')


if __name__ == '__main__':
    main()
