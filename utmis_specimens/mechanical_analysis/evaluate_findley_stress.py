from __future__ import print_function

from collections import namedtuple
import pickle
import os
import sys

import numpy as np

from case_hardening_toolbox.abaqus_functions.odb_io_functions import read_field_from_odb
from case_hardening_toolbox.abaqus_functions.odb_io_functions import write_field_to_odb
from case_hardening_toolbox.abaqus_functions.create_empty_odb import create_empty_odb

from transformation_fatigue.materials.gear_materials import SS2506
from transformation_fatigue.materials.gear_materials import SteelData

from transformation_fatigue.multiaxial_fatigue.multiaxial_fatigue_criteria import Findley


MechanicalData = namedtuple('MechanicalData', ['odb_file_name', 'step_name', 'frame_number'])


def perform_effective_stress_analysis(mechanical_data, effective_stress=Findley,
                                      element_set_name=None,
                                      instance_name=None, cpus=12, pickle_file=None, results_odb_file=None,
                                      results_odb_step_name=None, results_odb_frame_number=None, material=SS2506):
    """
    :param mechanical_data:             odb-file where mechanical simulations are found
    :param effective_stress:            The effective stress criterion to be used, default is Findley
    :param element_set_name             Element set where the effective stress should be evaluated. Default is None
                                        which gives the whole model
    :param instance_name                Name of the instance where data should be read. This parameter has to be
                                        specified if the odb consist of multiple instances. Default is None which apply
                                        to a model with one instance and an element set defined on assembly level
    :param cpus                         Number of cpus to run the effective fatigue stress evaluation, default is 8
    :param pickle_file                  Name of the pickle_file where the fatigue_data could be saved, default is
                                        None which does not save the data to a file
    :param results_odb_file             Odb file to write the results, default is None which does not write any data
                                        to a odb-file. If the file does not exist, an empty odb will be created using
                                        the heat_simulation_odb as base
    :param results_odb_step_name        Name of the step to write data, must be provided if results_odb_file is
                                        specified. An exception is thrown if it is not provided in such cases
    :param results_odb_frame_number     Frame number to write the data, default is None which creates a new frame at
                                        results_odb_step_name
    :param material                     Material to find material properties. Default is SS2506
    :return:                            A numpy array with the Fatigue stresses as the first column and
                                        Fatigue stress / critical fatigue stress as second column

    """
    if results_odb_file is not None and results_odb_step_name is None:
        raise ValueError("No results_odb_step_name is provided but results_odb_file is given")

    if results_odb_file is not None and not os.path.isfile(results_odb_file):
        create_empty_odb(results_odb_file, mechanical_data[0].odb_file_name)

    hardness_hv = read_field_from_odb('SDV_HARDNESS', mechanical_data[0].odb_file_name, set_name=element_set_name,
                                      instance_name=instance_name)

    stress_history = None
    for i, mechanical_data_set in enumerate(mechanical_data):
        stress = read_field_from_odb('S', mechanical_data_set.odb_file_name, step_name=mechanical_data_set.step_name,
                                     frame_number=mechanical_data_set.frame_number, set_name=element_set_name,
                                     instance_name=instance_name)
        if stress_history is None:
            stress_history = np.zeros((len(mechanical_data), stress.shape[0], 6))
        stress_history[i, :, :] = stress

    fatigue_data = np.zeros((stress_history.shape[1], 2))

    steel_data = SteelData(HV=hardness_hv)
    parameter = material.mean_stress_sensitivity(steel_data, effective_stress)
    critical_stress = material.critical_effective_stress(steel_data, effective_stress)
    fatigue_data[:, 0] = effective_stress.evaluation_function(stress_history, parameter, cpus=cpus)
    fatigue_data[:, 1] = fatigue_data[:, 0]/critical_stress

    if pickle_file:
        with open(pickle_file, 'w') as fatigue_pickle:
            pickle.dump(fatigue_data, fatigue_pickle)

    write_field_to_odb(field_data=fatigue_data[:, 0], field_id='SF', odb_file_name=results_odb_file,
                       step_name=results_odb_step_name, frame_number=results_odb_frame_number,
                       field_description=effective_stress.name, set_name=element_set_name,
                       instance_name=instance_name)
    frame_number = -1 if results_odb_frame_number is None else results_odb_frame_number
    write_field_to_odb(field_data=fatigue_data[:, 1], field_id='SFI', odb_file_name=results_odb_file,
                       step_name=results_odb_step_name, frame_number=frame_number,
                       field_description=(effective_stress.name + ' stress divided by critical '
                                          + effective_stress.name + ' stress'), set_name=element_set_name,
                       instance_name=instance_name)
    return fatigue_data


def main():
    k_450 = 0.3
    k_750 = 0.675
    b = (k_750 - k_450)/250
    a = k_750 - b*750
    SS2506.mean_stress_sensitivity_parameters = (a, b)
    specimen = sys.argv[-2]
    R = float(sys.argv[-1])
    specimen_loads = {'smooth': {-1.: [760 - 70, 760, 760 + 70], 0.: [424. - 26, 424, 424 + 26]},
                      'notched': {-1.: [439 - 20, 439, 439 + 20], 0.: [237. - 16, 237., 237 + 16]}}
    for load_amplitude in specimen_loads[specimen][R]:
        odb_file_directory = os.path.expanduser('~/utmis_specimens/' + specimen + '/mechanical_analysis/disp_control')
        sim_name = "snom=" + str(int(load_amplitude)) + "_R=" + str(int(R))
        odb_file_name = (odb_file_directory + "/utmis_" + specimen + '_' + sim_name + '.odb')

        results_odb_file = os.path.expanduser('~/utmis_specimens/' + specimen + '/findley_results_k='
                                              + str(k_750).replace('.', '_') + '.odb')

        mechanical_odb_data = [MechanicalData(odb_file_name=odb_file_name, step_name='3_max_load', frame_number=-1),
                               MechanicalData(odb_file_name=odb_file_name, step_name='3_min_load', frame_number=-1)]
        perform_effective_stress_analysis(mechanical_odb_data, element_set_name='FATIGUE_ELEMENTS',
                                          instance_name='SPECIMEN_PART_NEG', results_odb_file=results_odb_file,
                                          results_odb_step_name=sim_name)
        if R == -1.:
            perform_effective_stress_analysis(mechanical_odb_data, element_set_name='FATIGUE_ELEMENTS',
                                              instance_name='SPECIMEN_PART_POS', results_odb_file=results_odb_file,
                                              results_odb_step_name=sim_name, results_odb_frame_number=-1)


if __name__ == '__main__':
    main()
