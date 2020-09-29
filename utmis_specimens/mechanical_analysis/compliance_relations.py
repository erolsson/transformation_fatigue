import os

from mechanical_analysis_functions import write_mechanical_input_files, Simulation, Step, write_run_file
from transformation_fatigue.materials.materials import SS2506

def main():
    max_nominal_stress = 1000
    Wb = 2*5**2/6
    M = max_nominal_stress*Wb
    steps = [Step(name='loading', load=M)]
    simulations = [Simulation(name='loading', steps=steps, mode='force')]
    for specimen in ['smooth', 'notched']:
        geom_filename = os.path.expanduser('~/python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                           + specimen + '/utmis_' + specimen + '.inc')
        simulation_directory = os.path.expanduser('~/utmis_specimens/' + specimen + '/mechanical_analysis/')
        if not os.path.isdir(simulation_directory):
            os.makedirs(simulation_directory)
        job_names = write_mechanical_input_files(specimen, geom_filename, simulation_directory, simulations, SS2506)
        heat_treatment_data_file = os.path.expanduser('~/utmis_specimens/' + specimen + '/Toolbox_Cooling_utmis_'
                                                      + specimen + '.htd')
        write_run_file(job_names, heat_treatment_data_file, simulation_directory + '/run_compliance.sh')


if __name__ == '__main__':
    main()
