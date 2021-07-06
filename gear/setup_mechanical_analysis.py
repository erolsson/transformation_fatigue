import pathlib
import shutil
import sys

from transformation_fatigue.materials.materials import SS2506_02
from transformation_fatigue.transformation_subroutine.subroutine_info import subroutine_directory

from write_pulsator_model import create_pulsator_model, LoadStep


def create_input_files(simulation_directory):
    case_depth = float(sys.argv[-1])
    load_steps = 3
    R = 0.1  # noqa
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    number_of_teeth = 10
    input_files = {"dense":  model_directory / 'quarter_tooth_tilt2.inp',
                   "coarse": model_directory / 'coarse_tooth.inp'}
    load_amplitudes = [35.]
    simulation_names = []
    for load_amplitude in load_amplitudes:
        mean_load = (1 + R)/(1 - R)*load_amplitude
        steps = []
        for i in range(load_steps):
            steps.append(LoadStep(name="loading_" + str(i), load=(mean_load + load_amplitude)*1e3,
                                  output_time_interval=1.))
            steps.append(LoadStep(name="unloading_" + str(i), load=(mean_load - load_amplitude)*1e3,
                                  output_time_interval=1.))
        simulation_name = ("cd=" + str(case_depth).replace('.', '')
                           + "_Pamp=" + str(load_amplitude).replace('.', '_') + "kN")
        create_pulsator_model(simulation_name, input_files, simulation_directory, model_directory, number_of_teeth,
                              SS2506_02, steps)
        simulation_names.append(simulation_name)
    return simulation_names


def write_run_file(simulation_directory, simulation_names, pbs_job=True, cpus=12):
    file_lines = ["#!/bin/bash -i"]
    if pbs_job:
        file_lines.extend(['#PBS -V',
                           '#PBS -z',
                           '#PBS -l select=1:ncpus=' + str(cpus),
                           'cd $PBS_O_WORKDIR'])
    file_lines.append("shopt -s expand_aliases")
    file_lines.append("cd " + str(subroutine_directory))
    file_lines.append("./compile_gear_subroutine.sh")
    file_lines.append("cp gear_subroutine.o " + str(simulation_directory))
    file_lines.append("cd " + str(simulation_directory))
    for name in simulation_names:
        file_lines.append('abq2018 j=' + name + ' interactive cpus=' + str(cpus) + '  user=gear_subroutine.o')

    run_file_name = 'run_simulations_cd=' + sys.argv[-1].replace('.', '') + '.sh'
    with open(simulation_directory / run_file_name, 'w') as run_file:
        for line in file_lines:
            run_file.write(line + "\n")


def main():
    case_depth = float(sys.argv[-1])
    simulation_directory = (pathlib.Path.home() / "scania_gear_analysis" / "mechanical_analysis"
                            / ("cd=" + str(case_depth).replace('.', '')))
    if not simulation_directory.is_dir():
        simulation_directory.mkdir(parents=True)
    for tooth in ["coarse", "dense"]:
        heat_treatment_file = (pathlib.Path.home() / "scania_gear_analysis" / "odb_files" / "heat_treatment" /
                               ("cd=" + str(case_depth).replace('.', '') + '_' + tooth + '.htd'))
        shutil.copy(heat_treatment_file, simulation_directory)
    simulations = create_input_files(simulation_directory)
    write_run_file(simulation_directory, simulations)


if __name__ == '__main__':
    main()
