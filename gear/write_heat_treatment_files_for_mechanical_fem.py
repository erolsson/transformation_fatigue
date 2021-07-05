import pathlib

import numpy as np

from abaqus_python.abaqus_interface import ABQInterface

abq = ABQInterface('abq2018', output=True)
heat_sim_fields = ['AUSTENITE', 'CARBON', 'FERRITE', 'HARDNESS', 'LBAINITE', 'PEARLITE', 'Q_MARTENSITE',
                   'T_MARTENSITE', 'UBAINITE']


def main():
    for model_type in ['coarse', 'dense']:
        odb_directory = pathlib.Path.home() / "scania_gear_analysis" / "odb_files" / "heat_treatment"
        heat_treatment_odb = odb_directory / ('carbon_transfer_decarburization_' + model_type + '.odb')
        for cd in [0.5, 0.8, 1.1, 1.4]:
            step_name = 'results_cd=' + str(cd).replace('.', '')
            stress, _, element_labels = abq.read_data_from_odb('S', heat_treatment_odb, step_name=step_name,
                                                               get_position_numbers=True)
            num_points = stress.shape[0]
            num_elems = np.unique(element_labels).shape[0]
            num_gp = num_points//num_elems
            data = np.zeros((num_points, stress.shape[1] + 2 + len(heat_sim_fields)))
            gp = np.array(list(range(1, num_gp + 1, 1))*(num_points//num_gp))
            data[:, 0] = element_labels
            data[:, 1] = gp
            data[:, 2:2+stress.shape[1]] = stress
            for i, field_var in enumerate(heat_sim_fields):
                if field_var == 'HARDNESS':
                    field_var = 'HV'
                field = abq.read_data_from_odb(field_var, heat_treatment_odb, step_name=step_name)
                data[:, 2+stress.shape[1] + i] = field

            fmt = '%d, %d' + ', %f'*(data.shape[1] - 2)
            np.savetxt(odb_directory / ("cd=" + str(cd).replace('.', '') + '_' + model_type + '.htd'),
                       data, fmt=fmt, delimiter=',')


if __name__ == '__main__':
    main()
