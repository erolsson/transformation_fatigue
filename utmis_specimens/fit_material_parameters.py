import pathlib

from fat_eval.multiaxial_fatigue.criteria import haigh
from fat_eval.materials.fatigue_materials import Material
from abaqus_interface import ABQInterface

abq = ABQInterface('abq2018', output=True)


class SpecimenData:
    def __init__(self, pf_exp):
        self.stress_history = {}
        self.hardness = {}
        self.pf_exp = pf_exp


def evaluate_model(mechanical_data, effective_stress=haigh,
                   element_set_name=None,
                   instance_name=None, cpus=12, pickle_file=None, results_odb_file=None,
                   results_odb_step_name=None, results_odb_frame_number=None, material=SS2506):
    pass


def main():
    # specimen_loads = {'smooth': {-1.: [760 - 70, 760, 760 + 70], 0.: [424. - 26, 424, 424 + 26]},
    #                    'notched': {-1.: [439 - 20, 439, 439 + 20], 0.: [237. - 16, 237., 237 + 16]}}
    specimen_loads = {'smooth': {-1.: [], 0.: [424. - 26]},
                      'notched': {-1.: [], 0.: []}}
    data = {'smooth': {'-1': [], '0': []},
            'notched': {'-1': [], '0': []}}
    for specimen, load_data in specimen_loads.items():
        odb_file_directory = (pathlib.Path.home() / "utmis_specimens" / specimen / "mechanical_analysis_relaxed"
                              / "displacement_control_ra20")
        for r, loads in load_data.items():
            pf_exp = 0.25
            for load in loads:
                odb_filename = odb_file_directory / ("utmis_" + specimen + "_snom=" + str(int(load)) + "_R=" + str(int(r))
                                                     + ".odb")
                if 'elements' not in specimen_loads[specimen]:
                    stress, _, elements = abq.read_data_from_odb('SF', odb_filename, step_name='3_max_load',
                                                                 set_name='FATIGUE_ELEMENTS',
                                                                 instance_name='SPECIMEN_PART_NEG',
                                                                 get_position_numbers=True)
                    data[specimen]['elements'] = elements
                else:
                    stress = abq.read_data_from_odb('SF', odb_filename, step_name='3_max_load',
                                                    set_name='FATIGUE_ELEMENTS', instance_name='SPECIMEN_PART_NEG')
                stress_history = np.zeros((2, stress.shape[0], stress.shape[1]))
                stress_history[0, :, :] = stress
                specimen_data = SpecimenData(pf_exp)
                specimen_data.stress_history['neg'] =
                print(odb_file)
                pf_exp += 0.25


if __name__ == '__main__':
    main()
