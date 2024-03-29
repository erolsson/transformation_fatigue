import pathlib

from abaqus_python.abaqus_interface import ABQInterface, cylindrical_system_z, OdbInstance

from input_file_reader.input_file_reader import InputFileReader

from transformation_fatigue.materials.hardess_convertion_functions import HRC2HV

abq = ABQInterface('abq2018')
heat_sim_fields = ['AUSTENITE', 'CARBON', 'FERRITE', 'HARDNESS', 'LBAINITE', 'PEARLITE', 'Q_MARTENSITE',
                   'T_MARTENSITE', 'UBAINITE']


def main():
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    inp_file = model_directory / 'quarter_tooth_tilt2.inp'
    odb_directory = pathlib.Path.home() / "scania_gear_analysis" / "odb_files" / "heat_treatment"
    heat_sim_odb = odb_directory / 'carbon_transfer_decarburization.odb'
    dense_tooth = InputFileReader()
    dense_tooth.read_input_file(inp_file)
    abq.create_empty_odb_from_nodes_and_elements(heat_sim_odb, [OdbInstance('dense_tooth', dense_tooth)])

    heat_treatment_directory = (pathlib.Path.home() / "scania_gear_analysis" / "heat_simulation_dante_3"
                                / "carbon_transfer_decarburization")
    case_depths = [0.5, 0.8, 1.1, 1.4]
    for cd in case_depths:
        odb_file = heat_treatment_directory / ("CD" + str(cd).replace('.', '')) / "Toolbox_Cooling_VBC_gear.odb"
        step_name = 'results_cd=' + str(cd).replace('.', '')
        for field in heat_sim_fields:
            data = abq.read_data_from_odb(field_id='SDV_' + field, odb_file_name=odb_file)
            if field == 'HARDNESS':
                data = HRC2HV(data)
                field = 'HV'
            abq.write_data_to_odb(data, field, heat_sim_odb, step_name=step_name, frame_number=0)
        stress = abq.read_data_from_odb(field_id='S', odb_file_name=odb_file, coordinate_system=cylindrical_system_z)
        abq.write_data_to_odb(stress, 'S', heat_sim_odb, step_name=step_name, frame_number=0)


if __name__ == '__main__':
    main()
