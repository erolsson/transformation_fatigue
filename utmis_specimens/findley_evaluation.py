from __future__ import print_function

import os
import sys

import odbAccess

from case_hardening_toolbox.abaqus_functions.odb_io_functions import add_element_set


def get_list_from_set_file(filename):
    label_list = []
    with open(filename, 'r') as set_file:
        lines = set_file.readlines()
        for line in lines:
            label_list += line.split(',')
    return [int(label) for label in label_list]


def main():
    specimen = sys.argv[-2]
    R = int(sys.argv[-1])

    simulation_path = os.path.expanduser('~/utmis_specimens/')
    specimen_loads = {'smooth': {-1.: [737., 774., 820.], 0.: [425., 440.]},
                      'notched': {-1.: [427., 450.], 0.: [225., 240., 255.]}}

    for load in specimen_loads[specimen][R]:
        set_file = os.path.expanduser('~/python_projects/python_fatigue/fatigue_specimens/UTMIS/'
                                      'utmis_' + specimen + '/fatigue_volume_elements_' + specimen + '.inc')
        element_labels = get_list_from_set_file(set_file)
        odb_name = 'utmis_' + specimen + '_' + str(load).replace('.', '_') + '_R=' + str(R) + '.odb'
        odb_file = simulation_path + specimen + '/mechanical_analysis/' + odb_name
        add_element_set(odb_file, 'fatigue_volume_elements', element_labels, instance_name='SPECIMEN_PART_POS')
        add_element_set(odb_file, 'fatigue_volume_elements', element_labels, instance_name='SPECIMEN_PART_NEG')


if __name__ == '__main__':
    main()
