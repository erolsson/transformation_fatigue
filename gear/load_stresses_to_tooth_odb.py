import pathlib

from abaqus_python.abaqus_interface import ABQInterface, OdbInstance, cylindrical_system_z
from input_file_reader.input_file_reader import InputFileReader
from input_file_reader.input_file_functions import mirror_model


abq = ABQInterface('abq2018')


def create_tooth_odb(input_file, odb_file_name):
    positive_tooth = InputFileReader()
    positive_tooth.read_input_file(input_file)
    positive_tooth.renumber_nodes_and_elements()
    negative_tooth = mirror_model(positive_tooth, 'x')
    positive_instance = OdbInstance(name="positive_tooth", input_file_data=positive_tooth)
    negative_instance = OdbInstance(name="negative_tooth", input_file_data=negative_tooth)
    abq.create_empty_odb_from_nodes_and_elements(odb_file_name, [positive_instance, negative_instance])


def main():
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    inp_file = model_directory / 'quarter_tooth_tilt2.inp'
    tooth_odb_directory = pathlib.Path.home() / "scania_gear_analysis" / "odb_files" / "mechanical_analysis"
    tooth_odb_filename = tooth_odb_directory / "elastic_stresses.odb"
    tooth_odb_directory.mkdir(parents=True, exist_ok=True)
    create_tooth_odb(inp_file, tooth_odb_filename)

    mechanical_odb_filename = (pathlib.Path.home() / "scania_gear_analysis" / "mechanical_analysis" / "elastic"
                               / "pulsator_simulation.odb")
    for frame in abq.get_frames(mechanical_odb_filename):

        stress = abq.read_data_from_odb('S', odb_file_name=mechanical_odb_filename, frame_number=frame,
                                        coordinate_system=cylindrical_system_z, instance_name='EVAL_TOOTH_NEG')
        abq.write_data_to_odb(stress, 'S', tooth_odb_filename, 'loading', instance_name='negative_tooth')


if __name__ == '__main__':
    main()
