import os
import sys

import numpy as np

from mechanical_analysis_functions import write_mechanical_input_files, Simulation, Step, write_run_file
from transformation_fatigue.materials.materials import SS2506


def main():
    no_steps = 5
    specimen = sys.argv[-2]
    R = float(sys.argv[-1])
    specimen_loads = {'smooth': {-1.: [737., 774., 820.], 0.: [425., 440.]},
                      'notched': {-1.: [427., 450.], 0.: [225., 240., 255.]}}

    compliance_data = np.genfromtxt("compliance_utmis_" + specimen + ".csv", delimiter=",")
    stress_levels = compliance_data[:, 1]
    rotation = compliance_data[:, 2]
    simulations = []
    for load_ampiltude in specimen_loads[specimen][R]:
        mean_load = (1 + R)/(1 - R)*load_ampiltude
        max_load = mean_load + load_ampiltude
        mean_rotation = np.interp(mean_load, stress_levels, rotation)
        max_rotation = np.interp(max_load, stress_levels, rotation)

        rotation_amplitude = max_rotation - mean_rotation
        steps = []
        for step in range(1, no_steps + 1):
            steps.append(Step(str(step) + "_max_load", max_rotation, output_frequency=10.))
            steps.append(Step(str(step) + "_min_load", max_rotation - 2*rotation_amplitude, output_frequency=10.))
        simulations.append(Simulation("snom=" + str(int(load_ampiltude)) + "_R=" + str(int(R)), steps, 'displacement'))
    geom_filename = os.path.expanduser('~/python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                       + specimen + '/utmis_' + specimen + '.inc')
    simulation_directory = os.path.expanduser('~/utmis_specimens/' + specimen + '/mechanical_analysis/')
    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)
    job_names = write_mechanical_input_files(specimen, geom_filename, simulation_directory, simulations, SS2506)
    heat_treatment_data_file = os.path.expanduser('~/utmis_specimens/' + specimen + '/Toolbox_Cooling_utmis_'
                                                  + specimen + '.htd')
    write_run_file(job_names, heat_treatment_data_file, simulation_directory + '/run_R=' + str(int(R)) + '.sh')


if __name__ == '__main__':
    main()