import os
import sys

import numpy as np

from case_hardening_toolbox.abaqus_functions.odb_io_functions import read_field_from_odb
from case_hardening_toolbox.abaqus_functions.odb_io_functions import write_field_to_odb
from case_hardening_toolbox.abaqus_functions.create_empty_odb import create_empty_odb


def main():
    specimen_loads = {'smooth': {-1.: [760 - 70, 760, 760 + 70], 0.: [424. - 26, 424, 424 + 26]},
                      'notched': {-1.: [439 - 20, 439, 439 + 20], 0.: [237. - 16, 237., 237 + 16]}}
    specimen = sys.argv[-1]
    odb_file_directory = os.path.expanduser('~/utmis_specimens/' + specimen
                                            + '/mechanical_analysis/disp_control')
    results_odb_file = os.path.expanduser('~/utmis_specimens/' + specimen + '/stress_state.odb')
    sim_name = "snom=" + str(int(specimen_loads[specimen][-1][0])) + "_R=-1"
    create_empty_odb(results_odb_file, odb_file_directory + '/' + sim_name + '.odb')
    for R in [-1, 0]:
        for load_amplitude in specimen_loads[specimen][R]:
            sim_name = "snom=" + str(int(load_amplitude)) + "_R=" + str(int(R))
            odb_file_name = "/utmis_" + specimen + '_' + sim_name

            for instance_name in ['SPECIMEN_PART_NEG', 'SPECIMEN_PART_POS']:
                s1 = read_field_from_odb('S', odb_file_directory + '/' + odb_file_name + '.odb',
                                         step_name='3_min_load', frame_number=-1, instance_name=instance_name)

                s2 = read_field_from_odb('S', odb_file_directory + '/' + odb_file_name + '.odb',
                                         step_name='3_max_load', frame_number=-1, instance_name=instance_name)

                s_mean = (s1 + s2)/2
                s_amp = np.abs(s1 - s2)/2
                write_field_to_odb(field_data=s_mean, field_id='S_MEAN', odb_file_name=results_odb_file,
                                   step_name=sim_name, frame_number=None, field_description='mean stress',
                                   instance_name=instance_name)
                write_field_to_odb(field_data=s_amp, field_id='S_AMP', odb_file_name=results_odb_file,
                                   step_name=sim_name, frame_number=-1, field_description='amplitude stress',
                                   instance_name=instance_name)


if __name__ == '__main__':
    main()
