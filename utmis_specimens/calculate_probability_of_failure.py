from collections import namedtuple
import pathlib
import pickle

from abaqus_python.abaqus_interface import ABQInterface

SimulationData = namedtuple('SimulationData', ['stress', 'HV'])

abq = ABQInterface('abq2018', output=True)

specimen_loads = {'smooth': {-1.: [760 - 70, 760, 760 + 70], 0.: [424. - 26, 424, 424 + 26]},
                  'notched': {-1.: [439 - 20, 439, 439 + 20], 0.: [237. - 16, 237., 237 + 16]}}

ra = 15
k = 2.

odb_directory = pathlib.Path.home() / "utmis_specimens" / ("ra=" + str(int(ra)) + '_k=' + str(k).replace('.', '_'))


def load_findley_stress_states(pickle_file=None, write_pickle=False):
    if pickle_file and not write_pickle and pickle_file.is_file():
        with open(pickle_file, 'rb') as pickle_handle:
            return pickle.load(pickle_handle)

    data = {'smooth': {'-1': {'pos': {}, 'neg': {}}, '0': {'neg': {}}},
            'notched': {'-1': {'pos': {}, 'neg': {}}, '0': {'neg': {}}}}
    for specimen, r_data in specimen_loads.items():
        odb_filename = odb_directory / (specimen + '.odb')
        for r, load_levels in r_data.items():
            for load in load_levels:
                step_name = "snom=" + str(int(load)) + "_R=" + str(int(r))
                if 'elements' not in specimen_loads[specimen]:
                    stress, _, elements = abq.read_data_from_odb('SF', odb_filename, step_name=step_name,
                                                                 set_name='FATIGUE_ELEMENTS',
                                                                 instance_name='SPECIMEN_PART_NEG',
                                                                 get_position_numbers=True)
                    data[specimen]['elements'] = elements

                else:
                    stress = abq.read_data_from_odb('SF', odb_filename, step_name=step_name,
                                                    set_name='FATIGUE_ELEMENTS', instance_name='SPECIMEN_PART_NEG')
                hardness = abq.read_data_from_odb('HV', odb_filename, step_name=step_name,
                                                  set_name='FATIGUE_ELEMENTS', instance_name='SPECIMEN_PART_NEG')
                data[specimen][str(int(r))]['neg'][load] = SimulationData(stress, hardness)
                if r == -1:
                    stress = abq.read_data_from_odb('SF', odb_filename, step_name=step_name,
                                                    set_name='FATIGUE_ELEMENTS', instance_name='SPECIMEN_PART_POS')
                    hardness = abq.read_data_from_odb('HV', odb_filename, step_name=step_name,
                                                      set_name='FATIGUE_ELEMENTS', instance_name='SPECIMEN_PART_POS')
                    data[specimen][str(int(r))]['pos'][load] = SimulationData(stress, hardness)
    if write_pickle:
        with open(pickle_file, 'wb') as pickle_handle:
            pickle.dump(data, pickle_handle)
    return data


def main():
    fatigue_data = load_findley_stress_states(odb_directory / "data.pkl", True)
    pass


if __name__ == '__main__':
    main()
