import pathlib

import numpy as np
from scipy.optimize import fmin

from abaqus_interface import ABQInterface
from fat_eval.multiaxial_fatigue.evaluation import evaluate_effective_stress
from fat_eval.multiaxial_fatigue.criteria import haigh
from fat_eval.materials.fatigue_materials import Material
from fat_eval.weakest_link.weakest_link_evaluation import WeakestLinkEvaluator
from fat_eval.weakest_link.FEM_functions.node import Node
from fat_eval.weakest_link.FEM_functions.elements import element_types
from input_file_reader.input_file_reader import InputFileReader
from multiprocesser import multi_processer

abq = ABQInterface('abq2018', output=True)


class SpecimenData:
    def __init__(self, specimen, pf_exp):
        self.specimen = specimen
        self.stress_history = {}
        self.hardness = {}
        self.pf_exp = pf_exp
        self.elements = {}
        self.weakest_link_evaluators = {}

    def setup_weakest_link_evaluators(self):
        for part in self.stress_history:
            geom_filename = pathlib.Path.home() / ('python_projects/python_fatigue/fatigue_specimens/UTMIS/utmis_'
                                                   + self.specimen + '/utmis_' + self.specimen + '.inc')
            reader = InputFileReader()
            reader.read_input_file(geom_filename)
            nodes = dict(zip(reader.nodal_data[:, 0], reader.nodal_data[:, 1:]))
            model_elements = {}
            elements = []
            for element_type, element_data in reader.elements.items():
                model_elements.update(zip(element_data[:, 0], element_data[:, 1:]))
            for label in self.elements[part][::8]:
                element_nodes = [Node(coordinates=nodes[node], label=node) for node in model_elements[label]]
                elements.append(element_types['C3D8'](element_nodes))
            hv = self.hardness[part]
            hv = hv.reshape((hv.shape[0]//8, 8))
            self.weakest_link_evaluators[part] = WeakestLinkEvaluator(elements, hv, symmetry_factor=2)


def evaluate_model(specimen, material):
    ps = 1.
    for part in specimen.stress_history:
        s = evaluate_effective_stress(specimen.stress_history[part], material, criterion=haigh, cpus=12,
                                      hv=specimen.hardness[part])
        s = s.reshape(s.shape[0]//8, 8)
        ps *= (1 - specimen.weakest_link_evaluators[part].evaluate(s, material))
    return 1 - ps


def load_specimen_data(specimen, r, load, pf):
    odb_file_directory = (pathlib.Path.home() / "utmis_specimens" / specimen / "mechanical_analysis_relaxed"
                          / "displacement_control_ra20")
    odb_filename = odb_file_directory / ("utmis_" + specimen + "_snom=" + str(int(load)) + "_R="
                                         + str(int(r)) + ".odb")
    specimen_data = SpecimenData(specimen, pf)
    specimen_parts = ['neg']
    if r == -1:
        specimen_parts.append('pos')
    for part in specimen_parts:
        print("Starting reading data", odb_filename)
        hv, _, elements = abq.read_data_from_odb('SDV_HARDNESS', odb_filename, step_name='3_max_load',
                                                 set_name='FATIGUE_ELEMENTS',
                                                 instance_name='SPECIMEN_PART_' + part.upper(),
                                                 get_position_numbers=True)
        print("Hardness read")
        specimen_data.elements[part] = elements
        stress_history = np.zeros((2, hv.shape[0], 6))
        for i, step in enumerate(['max', 'min']):
            stress = abq.read_data_from_odb('S', odb_filename, step_name='3_' + step + '_load',
                                            set_name='FATIGUE_ELEMENTS',
                                            instance_name='SPECIMEN_PART_' + part.upper())
            stress_history[i, :, :] = stress

        specimen_data.stress_history[part] = stress_history
        specimen_data.hardness[part] = hv
    specimen_data.setup_weakest_link_evaluators()
    return specimen_data


def residual(par, specimens):
    mat = Material()
    mat.mean_stress_sensitivity_k = par[0]
    mat.sw1 = par[1]
    mat.sw2 = par[2]
    mat.b = par[3]
    pf_sim = []
    pf_exp = []
    for specimen in specimens:
        pf_sim.append(evaluate_model(specimen, mat))
        pf_exp.append(specimen.pf_exp)
    print(par)
    print(pf_sim)
    r = np.sum((np.array(pf_sim) - np.array(pf_exp))**2)
    print("Residual", r)
    return r


def main():
    specimen_loads = {'smooth': {-1.: [760 - 70, 760, 760 + 70], 0.: [424. - 26, 424, 424 + 26]},
                      'notched': {-1.: [439 - 20, 439, 439 + 20], 0.: [237. - 16, 237., 237 + 16]}}
    job_list = []
    for specimen, load_data in specimen_loads.items():
        for r, loads in load_data.items():
            pf_exp = 0.25
            for load in loads:
                job_list.append((load_specimen_data, [specimen, r, load, pf_exp], {}))
                pf_exp += 0.25
    print("Stating collecting data")
    specimen_data = multi_processer(job_list, cpus=1, timeout=3600)
    fmin(residual, [1000, 1000, 0, 6e6], (specimen_data,))


if __name__ == '__main__':
    main()
