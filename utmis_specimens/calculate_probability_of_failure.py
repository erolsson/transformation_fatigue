from collections import namedtuple
import pathlib
import pickle

import numpy as np

from scipy.optimize import fmin

from abaqus_python.abaqus_interface import ABQInterface

from fat_eval.FEM_functions.node import Node
from fat_eval.FEM_functions.elements import element_types
from fat_eval.weakest_link.weakest_link_evaluation import WeakestLinkEvaluator
from fat_eval.materials.fatigue_materials import Material

from input_file_reader.input_file_reader import InputFileReader

from multiprocesser.multiprocesser import multi_processer

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


def setup_weakest_link_calculators(fatigue_data):
    evaluators = {'smooth': {'-1': {'pos': {}, 'neg': {}}, '0': {'neg': {}}},
                  'notched': {'-1': {'pos': {}, 'neg': {}}, '0': {'neg': {}}}}
    for specimen, r_data in fatigue_data.items():
        geom_filename = pathlib.Path.home() / ('python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                               + specimen + '/utmis_' + specimen + '.inc')
        reader = InputFileReader()
        reader.read_input_file(geom_filename)
        nodes = dict(zip(reader.nodal_data[:, 0], reader.nodal_data[:, 1:]))
        model_elements = {}
        elements = []
        for element_type, element_data in reader.elements.items():
            model_elements.update(zip(element_data[:, 0], element_data[:, 1:]))
        for label in fatigue_data[specimen]['elements'][::8]:
            element_nodes = [Node(coordinates=nodes[node], label=node) for node in model_elements[label]]
            elements.append(element_types['C3D8'](element_nodes))

        for r in ['-1', '0']:
            for instance, load_levels in r_data[r].items():
                for load in load_levels:
                    hv = fatigue_data[specimen][r][instance][load].HV
                    hv = hv.reshape(hv.shape[0]//8, 8)
                    evaluators[specimen][r][instance][load] = WeakestLinkEvaluator(elements, hv, symmetry_factor=2)
    return evaluators


def calculate_probability_of_failure(specimen, r, load, material, fatigue_data, evaluators):
    ps = 1.
    instances = ['neg']
    if r == -1:
        instances.append('pos')
    for instance in instances:
        evaluator = evaluators[specimen][str(r)][instance][load]
        stress = fatigue_data[specimen][str(r)][instance][load].stress
        stress = stress.reshape(stress.shape[0]//8, 8)
        ps *= (1-evaluator.evaluate(stress, material))
    return 1-ps


def evaluate_probabilities_of_failure(evaluators, fatigue_data, material):
    job_list = []
    for specimen, r_data in fatigue_data.items():
        for r in ['-1', '0']:
            for load in r_data[r]['neg']:
                job_list.append((calculate_probability_of_failure,
                                [specimen, r, load, material, fatigue_data, evaluators], {}))

    return multi_processer(job_list, delay=0, timeout=3600)


def fit_material_parameters(evaluators, fatigue_data):
    def residual(par):
        mat = Material()
        mat.sw1 = par[0]
        mat.sw2 = par[1]
        mat.b = par[2]
        pf_sim = evaluate_probabilities_of_failure(evaluators, fatigue_data, mat)
        pf_exp = np.array([0.5 - 0.341, 0.5, 0.5 + 0.341]*4)
        r = np.sum((pf_sim - pf_exp)**2)
        print(mat.sw1, mat.sw2, mat.b, r)
        print(pf_sim)
        return r
    print(fmin(residual, [1000, 0, 6e6]))


def main():
    fatigue_data = load_findley_stress_states(odb_directory / "data.pkl")
    evaluators = setup_weakest_link_calculators(fatigue_data)
    test_material = Material()
    test_material.b = 6e6
    test_material.sw1 = 1200
    fit_material_parameters(evaluators, fatigue_data)


if __name__ == '__main__':
    main()
