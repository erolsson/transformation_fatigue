import os
import sys

import numpy as np

from mechanical_analysis_functions import write_mechanical_input_files, Simulation, Step, write_run_file
from transformation_fatigue.materials.materials import SS2506 as SS2506


def main():
    wb = 2*25/6
    no_steps = 3
    specimen = sys.argv[-2]
    R = float(sys.argv[-1])
    # specimen_loads = {'smooth': {-1.: [760 - 70, 760, 760 + 70], 0.: [424. - 26, 424, 424 + 26]},
    #                   'notched': {-1.: [439 - 20, 439, 439 + 20], 0.: [237. - 16, 237., 237 + 16]}}
    specimen_loads = {'smooth': {-1.: [760], 0.: [424]}, 'notched': {-1.: [439], 0.: [237]}}
    heat_treatment_simulation = 't=9min_90C_decarburization'
    simulations = []
    compliance_data = np.genfromtxt('compliance_utmis_' + specimen + '.csv', delimiter=',')
    for load_amplitude in specimen_loads[specimen][R]:
        amplitude_rot = np.interp(load_amplitude, compliance_data[:, 1], compliance_data[:, 2])
        mean_rot = 0
        if R == 0.:
            mean_rot = amplitude_rot
        steps = []
        for step in range(1, no_steps + 1):
            steps.append(Step(str(step) + "_max_load", mean_rot + amplitude_rot, output_frequency=10.))
            steps.append(Step(str(step) + "_min_load", mean_rot - amplitude_rot, output_frequency=10.))
        simulations.append(Simulation("snom=" + str(int(load_amplitude)) + "_R=" + str(int(R)), steps, 'displacement'))
    geom_filename = os.path.expanduser('~/python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                       + specimen + '/utmis_' + specimen + '.inc')
    simulation_directory = os.path.expanduser('~/utmis_specimens/' + specimen
                                              + '/mechanical_analysis/disp_control/')
    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)
    job_names = write_mechanical_input_files(specimen, geom_filename, simulation_directory, simulations, SS2506)
    heat_treatment_data_file = os.path.expanduser('~/utmis_specimens/' + specimen + '/' + heat_treatment_simulation +
                                                  '/Toolbox_Cooling_utmis_' + specimen + '.htd')
    write_run_file(job_names, heat_treatment_data_file, simulation_directory + '/run_R=' + str(int(R)) + '.sh')


if __name__ == '__main__':
    main()
