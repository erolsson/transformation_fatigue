import pathlib

import numpy as np

from input_file_reader.input_file_functions import mirror_model
from input_file_reader.input_file_reader import InputFileReader
from transformation_fatigue.materials.materials import SS2506, SS2506Elastic


class GearTooth:
    def __init__(self, tooth_number, rotation, part_names, position=(0., 0., 0.)):
        self.rotation = rotation
        self.part_names = part_names
        self.pos = position
        self.instance_names = [part + '_' + str(tooth_number) for part in self.part_names]

    def write_input(self):
        lines = []
        for instance_name, part_name in zip(self.instance_names, self.part_names):
            lines += ['\t*Instance, name=' + instance_name + ', part=' + part_name,
                      '\t\t' + str(self.pos[0]) + ', ' + str(self.pos[1]) + ', ' + str(self.pos[2]),
                      '\t\t' + str(self.pos[0]) + ', ' + str(self.pos[1]) + ', ' + str(self.pos[2]) + ', ' +
                      str(self.pos[0]) + ', ' + str(self.pos[1]) + ', ' + str(self.pos[2] + 1.0) + ', ' +
                      str(self.rotation),
                      '\t*End Instance']
        return lines


class LoadStep:
    def __init__(self, name, load, time=1., initial_increment=0.01, output_time_interval=None):
        self.name = name
        self.load = load
        self.duration = time
        self.initial_increment = initial_increment
        self.output_time_interval = output_time_interval

    def write_input(self):
        lines = ['*step, name=' + self.name + ', nlgeom=Yes',
                 '\t*Static',
                 '\t\t' + str(self.initial_increment) + ', ' + str(self.duration) + ', 1e-12, 1.',
                 '\t*CLoad',
                 '\t\tjaw_ref_node, 1, -0.5',
                 '\t*CLoad',
                 '\t\tjaw_ref_node, 2, ' + str(-self.load/2)]
        if self.output_time_interval is None:
            lines.append('\t*Output, field')
        else:
            lines.append('\t*Output, field, time interval=' + str(self.output_time_interval))
        lines.append('\t\t*Element Output, directions=no')
        lines.append('\t\t\tS, SDV')
        lines.append('\t\t*Node Output')
        lines.append('\t\t\tU')
        lines.append('*End step')
        return lines


def write_tooth_part(name, inc_file, set_file):
    lines = ['*Part, name=' + name,
             '\t*Include, Input=' + inc_file,
             '\t*Include, Input=' + set_file,
             '\t*Solid Section, elset=all_elements, material=SS2506',
             '\t\t1.0',
             '*End Part']
    return lines


def create_jaw_sets(jaw_input_file_reader):
    jaw_nodes = jaw_input_file_reader.nodal_data
    jaw_elements = jaw_input_file_reader.elements
    y = np.unique(jaw_nodes[:, 2])
    x_min = np.unique(jaw_nodes[:, 1])[0]
    y_min, y_max = y[0], y[-1]
    node_sets = {'x_min_nodes': jaw_nodes[jaw_nodes[:, 1] == x_min, 0],
                 'y_min_nodes': jaw_nodes[jaw_nodes[:, 2] == y_min, 0],
                 'y_max_nodes': jaw_nodes[jaw_nodes[:, 2] == y_max, 0],
                 'z0_nodes': jaw_nodes[jaw_nodes[:, 3] == 0.0, 0]}

    y_min_elements = []
    y_max_elements = []
    x_min_elements = []

    for element_data in jaw_elements.values():
        for e in element_data:
            for n_label in e[1:]:
                for element_list, node_label_set in zip([x_min_elements, y_min_elements, y_max_elements],
                                                        [node_sets['x_min_nodes'], node_sets['y_min_nodes'],
                                                         node_sets['y_max_nodes']]):
                    if n_label in node_label_set:
                        element_list.append(e[0])

    element_sets = {'x_min_elements': np.unique(x_min_elements),
                    'y_min_elements': np.unique(y_min_elements),
                    'y_max_elements': np.unique(y_max_elements)}

    for set_name, nodes in node_sets.items():
        jaw_input_file_reader.create_node_set(set_name, nodes)

    for set_name, elements in element_sets.items():
        jaw_input_file_reader.create_element_set(set_name, elements)


def create_pulsator_model(job_name, input_files, simulation_directory, model_directory, number_of_teeth, material,
                          load_steps):
    for name, input_file in input_files.items():
        positive_tooth = InputFileReader()
        positive_tooth.read_input_file(input_file)
        positive_tooth.renumber_nodes_and_elements()
        positive_tooth.write_geom_include_file(simulation_directory / (name + "_tooth_pos.inc"))
        positive_tooth.write_sets_file(simulation_directory / (name + "_tooth_sets.inc"),
                                       surfaces_from_element_sets=['exposed', 'x0', 'x1'])
        negative_tooth = mirror_model(positive_tooth, 'x')
        negative_tooth.write_geom_include_file(simulation_directory / (name + "_tooth_neg.inc"))

    reader = InputFileReader()
    reader.read_input_file(model_directory / "pulsator_jaw.inp")
    z0 = np.min(np.abs(np.unique(reader.nodal_data[:, 3])))
    reader.nodal_data[:, 3] += z0
    jaw_nodes = reader.nodal_data[reader.nodal_data[:, 3] >= 0., :]
    reader.remove_nodes(jaw_nodes[:, 0])
    # Swap x and y
    temp = jaw_nodes[:, 1].copy()
    jaw_nodes[:, 1] = jaw_nodes[:, 2]
    jaw_nodes[:, 2] = temp
    jaw_nodes[:, 2] *= -1
    reader.nodal_data = jaw_nodes
    reader.renumber_nodes_and_elements()
    create_jaw_sets(reader)
    reader.write_geom_include_file(simulation_directory / "pulsator_jaw_geom.inc")
    reader.write_sets_file(simulation_directory / "pulsator_jaw_sets.inc",
                           surfaces_from_element_sets=['x_min', 'y_min', 'y_max'])

    teeth = []
    part_names = [('coarse_tooth_pos', 'coarse_tooth_neg') for _ in range(number_of_teeth)]
    part_names[0] = ('coarse_tooth_pos', 'dense_tooth_neg')
    part_names[1] = ('dense_tooth_pos', 'dense_tooth_neg')
    part_names[2] = ('dense_tooth_pos', 'coarse_tooth_neg')
    for i in range(number_of_teeth):
        teeth.append(GearTooth(tooth_number=i,
                               rotation=(180/number_of_teeth)*i - 90. + (180/2/number_of_teeth),
                               part_names=part_names[i]))

    # Tooth number 1 is the interesting tooth for fatigue, give it a denser mesh and a different name
    file_lines = ['*Heading',
                  '\tModel of a pulsator test of a planetary gear']
    for mesh in ['coarse', 'dense']:
        for sign in ['pos', 'neg']:
            file_lines += write_tooth_part(name=mesh + '_tooth_' + sign, inc_file=mesh + '_tooth_' + sign + '.inc',
                                           set_file=mesh + '_tooth_sets.inc')

    file_lines.append('*Part, name=pulsator_jaw_part')
    file_lines.append('\t*Include, input=pulsator_jaw_geom.inc')
    file_lines.append('\t*Include, input=pulsator_jaw_sets.inc')
    file_lines.append('\t*Solid Section, elset=all_elements, material=SS2506Elastic')
    file_lines.append('\t\t1.0')
    file_lines.append('*End Part')

    file_lines.append('**')
    file_lines.extend(material.material_input_file_string())
    file_lines.extend(SS2506Elastic.material_input_file_string())
    file_lines.append('*Assembly, name=pulsator_model')

    for tooth in teeth:
        file_lines += tooth.write_input()

    file_lines.append('\t*Instance, name=pulsator_jaw, part=pulsator_jaw_part')
    file_lines.append('\t*End Instance')

    # Writing the tie constraints at the mid lines of the teeth
    for i in range(number_of_teeth):
        file_lines.append('\t*Tie, name=tie_mid_tooth' + str(i))
        file_lines.append('\t\t' + teeth[i].instance_names[0] + '.x0_surface, '
                          + teeth[i].instance_names[1] + '.x0_surface')

    # Writing tie constraints between the teeth
    for i in range(1, number_of_teeth):
        file_lines.append('\t*Tie, name=tie_inter_teeth_' + str(i-1) + '_' + str(i))
        file_lines.append('\t\t' + teeth[i-1].instance_names[1] + '.x1_surface, '
                          + teeth[i].instance_names[0] + '.x1_surface')

    # Adding a kinematic coupling for the pulsator jaw
    file_lines.append('\t*Node, nset=jaw_ref_node')
    file_lines.append('\t\t999999, ' + str(np.min(jaw_nodes[:, 1])) + ',' + str(np.max(jaw_nodes[:, 2])) + ', 0.0')
    file_lines.append('\t*Nset, nset=pinned, instance=' + teeth[0].instance_names[0])
    file_lines.append('\t\t21')
    file_lines.append('\t*Coupling, Constraint name=jaw_load_coupling, '
                      'ref node=jaw_ref_node, surface=Pulsator_jaw.y_max_surface')
    file_lines.append('\t\t*Kinematic')
    file_lines.append('*End Assembly')

    # Creating the contact between the pulsator jaw and the eval tooth in the vertical direction
    file_lines.append('*Surface interaction, name=frictionless_contact')
    file_lines.append('*Contact pair, interaction=frictionless_contact')
    file_lines.append('\tdense_tooth_neg_1.exposed_surface, Pulsator_jaw.y_min_surface')
    file_lines.append('*Contact pair, interaction=frictionless_contact')
    file_lines.append('\tdense_tooth_pos_2.exposed_surface, Pulsator_jaw.x_min_surface')

    for tooth in teeth:
        file_lines.append('*Boundary')
        file_lines.append('\t' + tooth.instance_names[0] + '.z0_nodes, 3, 3')
        file_lines.append('*Boundary')
        file_lines.append('\t' + tooth.instance_names[1] + '.z0_nodes, 3, 3')

    file_lines.append('*Boundary')
    file_lines.append('\tpulsator_jaw.z0_nodes, 3, 3')

    file_lines.append('*Boundary')
    file_lines.append('\tpinned, 1, 1')

    file_lines.append('*Boundary')
    file_lines.append('\t' + teeth[0].instance_names[0] + '.x1_nodes, 2, 2')
    file_lines.append('*Boundary')
    file_lines.append('\t' + teeth[-1].instance_names[1] + '.x1_nodes, 2, 2')

    file_lines.append('*Boundary')
    file_lines.append('\tjaw_ref_node, 3, 6')
    file_lines.append('*Initial Conditions, type=Solution, user')
    file_lines.append('*Initial Conditions, type=Stress, user')
    file_lines.append('*Initial conditions, type=temperature')
    for tooth in teeth:
        file_lines.append('\t' + tooth.instance_names[0] + '.ALL_NODES, 20')
        file_lines.append('\t' + tooth.instance_names[1] + '.ALL_NODES, 20')

    initiate_contact_step = LoadStep('Initiate_contact', initial_increment=1e-4, load=1.)

    initiate_contact_lines = initiate_contact_step.write_input()
    initiate_contact_lines.insert(3, '\t*Controls, reset')
    initiate_contact_lines.insert(4, '\t*Controls, parameters=line search')
    initiate_contact_lines.insert(5, '\t\t5, , , , ')
    initiate_contact_lines.insert(6, '\t*Contact Controls, Stabilize')
    initiate_contact_lines.insert(7, '\t*Contact Interference, shrink')
    initiate_contact_lines.insert(8, '\t\tdense_tooth_neg_1.exposed_surface, Pulsator_jaw.y_min_surface')
    initiate_contact_lines.insert(9, '\t*Contact Interference, shrink')
    initiate_contact_lines.insert(10, '\t\tdense_tooth_pos_2.exposed_surface, Pulsator_jaw.x_min_surface')

    file_lines += initiate_contact_lines
    for load_step in load_steps:
        file_lines += load_step.write_input()

    with open(simulation_directory / (job_name + '.inp'), 'w') as input_file:
        for line in file_lines:
            input_file.write(line + '\n')
        input_file.write("**EOF")


def main():
    max_load = 130.
    simulation_directory = pathlib.Path.home() / "scania_gear_analysis" / "mechanical_analysis" / "cd=05"
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    number_of_teeth = 10
    input_files = {"dense":  model_directory / 'quarter_tooth_tilt2.inp',
                   "coarse": model_directory / 'coarse_tooth.inp'}

    load_steps = [LoadStep('loading', max_load*1e3, output_time_interval=0.02)]
    if not simulation_directory.is_dir():
        simulation_directory.mkdir(parents=True)
    create_pulsator_model("pulsator_simulation", input_files, simulation_directory, model_directory, number_of_teeth,
                          SS2506, load_steps)


if __name__ == '__main__':
    main()
