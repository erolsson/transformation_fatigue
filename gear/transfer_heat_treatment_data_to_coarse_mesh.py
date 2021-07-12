from collections import namedtuple
import pathlib

import numpy as np

from abaqus_interface.abaqus_interface import ABQInterface, OdbInstance
from input_file_reader.input_file_reader import InputFileReader

from FEM_functions.elements import C3D8
from create_heat_treatment_odb import heat_sim_fields

abq = ABQInterface('abq2018')

Node = namedtuple('Node', ['coordinates', 'label'])


def main():
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    coarse_input_file = model_directory / 'coarse_tooth.inp'
    coarse_tooth = InputFileReader()
    coarse_tooth.read_input_file(coarse_input_file)
    coarse_tooth.renumber_nodes_and_elements()
    odb_directory = pathlib.Path.home() / "scania_gear_analysis" / "odb_files" / "heat_treatment"
    dense_tooth_odb = odb_directory / 'carbon_transfer_decarburization_dense.odb'
    coarse_tooth_odb = odb_directory / 'carbon_transfer_decarburization_coarse.odb'
    abq.create_empty_odb_from_nodes_and_elements(coarse_tooth_odb, [OdbInstance('dense_tooth', coarse_tooth)])

    nodes = dict(zip(coarse_tooth.nodal_data[:, 0], coarse_tooth.nodal_data[:, 1:]))
    gauss_points = []
    for element_list in coarse_tooth.elements.values():
        for element in element_list:
            element_nodes = [Node(coordinates=nodes[node], label=node) for node in element[1:]]
            e = C3D8(element_nodes)
            for i in range(8):
                gauss_points.append(e.gauss_point_coordinates(i))
    gauss_points = np.array(gauss_points)

    for cd in [0.5, 0.8, 1.1, 1.4]:
        step_name = 'results_cd=' + str(cd).replace('.', '')
        for field_name in heat_sim_fields:
            if field_name == 'HARDNESS':
                field_name = 'HV'
            field = abq.get_data_from_path(gauss_points, dense_tooth_odb, variable=field_name, step_name=step_name,
                                           output_position='INTEGRATION_POINT')
            abq.write_data_to_odb(field, field_name, coarse_tooth_odb, step_name, frame_number=0)
        stress = abq.get_tensor_from_path(odb_file_name=dense_tooth_odb, path_points=gauss_points, field_id='S',
                                          step_name=step_name)
        abq.write_data_to_odb(stress, 'S', coarse_tooth_odb, step_name, frame_number=0)


if __name__ == '__main__':
    main()
