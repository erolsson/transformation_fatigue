import pathlib

from mechanical_analysis_functions import write_mechanical_input_files, Simulation, Step, write_run_file
from transformation_fatigue.materials.materials import SS2506_02 as SS2506


def main():
    max_nominal_stress = 1000
    Wb = 2*5**2/6
    M = max_nominal_stress*Wb
    steps_smooth = [Step(name='loading', load=M, max_inc=0.01)]
    steps_notched = [Step(name='loading', load=M/10*6, max_inc=0.01)]
    simulations = {"smooth": [Simulation(name='loading', steps=steps_smooth, mode='force')],
                   "notched": [Simulation(name='loading', steps=steps_notched, mode='force')]}
    for specimen in ['smooth', 'notched']:
        geom_filename = (pathlib.Path.home() / 'python_projects/python_fatigue/fatigue_specimens/UTMIS'
                         / ('utmis_' + specimen) / ('utmis_' + specimen + '.inc'))

        simulation_directory = (pathlib.Path.home() / 'utmis_specimens' / specimen / 'mechanical_analysis_relaxed_02')
        if not simulation_directory.is_dir():
            simulation_directory.mkdir(parents=True)
        job_names = write_mechanical_input_files(specimen, geom_filename, simulation_directory, simulations[specimen],
                                                 SS2506)
        heat_treatment_simulation = 't=9min_75C_decarburization'
        heat_treatment_data_file = (pathlib.Path.home() / 'utmis_specimens/' / specimen / heat_treatment_simulation
                                    / ('Toolbox_Cooling_utmis_' + specimen + '.htd'))
        write_run_file(job_names, heat_treatment_data_file, simulation_directory / 'run_compliance.sh')


if __name__ == '__main__':
    main()
