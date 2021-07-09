import os
import pathlib
import sys

import numpy as np

from mechanical_analysis_functions import write_mechanical_input_files, Simulation, Step, write_run_file
from transformation_fatigue.materials.materials import SS2506_02 as SS2506


def main():
    wb = 2*25/6
    no_steps = 3
    specimen = sys.argv[-2]
    R = float(sys.argv[-1])
    specimen_loads = {'smooth': {-1.: [760 - 70, 760, 760 + 70], 0.: [424. - 26, 424, 424 + 26]},
                      'notched': {-1.: [439 - 20, 439, 439 + 20], 0.: [237. - 16, 237., 237 + 16]}}
    # specimen_loads = {'smooth': {-1.: [760], 0.: [424]}, 'notched': {-1.: [439], 0.: [237]}}
    heat_treatment_simulation = 't=9min_75C_decarburization'
    simulations = []
    compliance_data = np.genfromtxt('compliance_utmis_' + specimen + '.csv', delimiter=',')
    loading = 'displacement'
    for load_amplitude in specimen_loads[specimen][R]:
        if loading == 'displacement':
            amplitude = np.interp(load_amplitude, compliance_data[:, 1], compliance_data[:, 2])
        else:
            amplitude = load_amplitude*wb
        mean = 0
        if R == 0.:
            mean = amplitude
        steps = []
        for step in range(1, no_steps + 1):
            steps.append(Step(str(step) + "_max_load", mean + amplitude, output_frequency=10.))
            steps.append(Step(str(step) + "_min_load", mean - amplitude, output_frequency=10.))
        simulations.append(Simulation("snom=" + str(int(load_amplitude)) + "_R=" + str(int(R)), steps, 'displacement'))
    geom_filename = pathlib.Path.home() / ('python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                           + specimen + '/utmis_' + specimen + '.inc')
    simulation_directory = pathlib.Path.home() / ('utmis_specimens/' + specimen
                                                  + '/mechanical_analysis_relaxed/' + loading + '_control_ra20/')
    if not os.path.isdir(simulation_directory):
        os.makedirs(simulation_directory)
    job_names = write_mechanical_input_files(specimen, geom_filename, simulation_directory, simulations, SS2506)
    heat_treatment_data_file = pathlib.Path.home() / ('utmis_specimens/' + specimen + '/' + heat_treatment_simulation +
                                                      '/Toolbox_Cooling_utmis_' + specimen + '.htd')
    write_run_file(job_names, heat_treatment_data_file, simulation_directory / ('run_R=' + str(int(R)) + '.sh'))


if __name__ == '__main__':
    main()
