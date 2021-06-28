import pathlib

import numpy as np

from input_file_reader.input_file_functions import mirror_model
from input_file_reader.input_file_reader import InputFileReader


class GearTooth:
    def __init__(self, instance_name, rotation, part_names, position=(0., 0., 0.)):
        self.instance_name = instance_name
        self.rotation = rotation
        self.part_names = part_names
        self.pos = position
        self.sign_names = ['pos', 'neg']

    def write_input(self):
        lines = []
        for sign, part_name in zip(self.sign_names, self.part_names):
            lines += ['\t*Instance, name=' + self.instance_name + '_' + str(sign) + ', part=' + part_name,
                      '\t\t' + str(self.pos[0]) + ', ' + str(self.pos[1]) + ', ' + str(self.pos[2]),
                      '\t\t' + str(self.pos[0]) + ', ' + str(self.pos[1]) + ', ' + str(self.pos[2]) + ', ' +
                      str(self.pos[0]) + ', ' + str(self.pos[1]) + ', ' + str(self.pos[2] + 1.0) + ', ' +
                      str(self.rotation),
                      '\t*End Instance']
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


def write_load_step(step_name, force=None, step_time=1., initial_inc=0.01, output_time_interval=None):
    lines = ['*step, name=' + step_name + ', nlgeom=Yes',
             '\t*Static',
             '\t\t' + str(initial_inc) + ', ' + str(step_time) + ', 1e-12, 1.',
             '\t*CLoad',
             '\t\tjaw_ref_node, 1, -0.5']

    if force:
        lines.append('\t*CLoad')
        lines.append('\t\tjaw_ref_node, 2, ' + str(-force/2))
    if output_time_interval is None:
        lines.append('\t*Output, field')
    else:
        lines.append('\t*Output, field, time interval=' + str(output_time_interval))
    lines.append('\t\t*Element Output')
    lines.append('\t\t\tS')
    lines.append('\t\t*Node Output')
    lines.append('\t\t\tU')
    lines.append('*End step')
    return lines


def main():
    max_load = 130.
    simulation_directory = pathlib.Path.home() / "scania_gear_analysis" / "mechanical_analysis" / "elastic"
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    number_of_teeth = 10
    input_files = {"dense":  model_directory / 'quarter_tooth_tilt2.inp',
                   "coarse": model_directory / 'coarse_tooth.inp'}

    if not simulation_directory.is_dir():
        simulation_directory.mkdir(parents=True)
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
    for i in range(number_of_teeth):
        teeth.append(GearTooth(instance_name='tooth' + str(i),
                               rotation=(180/number_of_teeth)*i - 90. + (180/2/number_of_teeth),
                               part_names=['coarse_tooth_pos', 'coarse_tooth_neg']))

    # Tooth number 1 is the interesting tooth for fatigue, give it a denser mesh and a different name
    teeth[1].instance_name = 'eval_tooth'
    teeth[0].part_names = ['coarse_tooth_pos', 'dense_tooth_neg']
    teeth[1].part_names = ['dense_tooth_pos', 'dense_tooth_neg']
    teeth[2].part_names = ['dense_tooth_pos', 'coarse_tooth_neg']

    file_lines = ['*Heading',
                  '\tModel of a pulsator test of a planetary gear']
    for mesh in ['coarse', 'dense']:
        for sign in ['pos', 'neg']:
            file_lines += write_tooth_part(name=mesh + '_tooth_' + sign, inc_file=mesh + '_tooth_' + sign + '.inc',
                                           set_file=mesh + '_tooth_sets.inc')

    file_lines.append('*Part, name=pulsator_jaw_part')
    file_lines.append('\t*Include, input=pulsator_jaw_geom.inc')
    file_lines.append('\t*Include, input=pulsator_jaw_sets.inc')
    file_lines.append('\t*Solid Section, elset=all_elements, material=SS2506')
    file_lines.append('\t\t1.0')
    file_lines.append('*End Part')

    file_lines.append('**')
    file_lines.append('*Material, name=SS2506')
    file_lines.append('\t*Elastic')
    file_lines.append('\t\t200E3, 0.3')
    file_lines.append('*Assembly, name=pulsator_model')

    for tooth in teeth:
        file_lines += tooth.write_input()

    file_lines.append('\t*Instance, name=pulsator_jaw, part=pulsator_jaw_part')
    file_lines.append('\t*End Instance')

    # Writing the tie constraints at the mid lines of the teeth
    for i in range(number_of_teeth):
        file_lines.append('\t*Tie, name=tie_mid_tooth' + str(i))
        file_lines.append('\t\t' + teeth[i].instance_name +
                          '_pos.x0_surface, ' + teeth[i].instance_name + '_neg.x0_surface')

    # Writing tie constraints between the teeth
    for i in range(1, number_of_teeth):
        file_lines.append('\t*Tie, name=tie_inter_teeth_' + str(i-1) + '_' + str(i))
        file_lines.append('\t\t' + teeth[i-1].instance_name +
                          '_neg.x1_surface, ' + teeth[i].instance_name + '_pos.x1_surface')

    # Adding a kinematic coupling for the pulsator jaw
    file_lines.append('\t*Node, nset=jaw_ref_node')
    file_lines.append('\t\t999999, ' + str(np.min(jaw_nodes[:, 1])) + ',' + str(np.max(jaw_nodes[:, 2])) + ', 0.0')
    file_lines.append('\t*Nset, nset=pinned, instance=' + teeth[0].instance_name + '_pos')
    file_lines.append('\t\t21')
    file_lines.append('\t*Coupling, Constraint name=jaw_load_coupling, '
                      'ref node=jaw_ref_node, surface=Pulsator_jaw.y_max_surface')
    file_lines.append('\t\t*Kinematic')
    file_lines.append('*End Assembly')

    # Creating the contact between the pulsator jaw and the eval tooth in the vertical direction
    file_lines.append('*Surface interaction, name=frictionless_contact')
    file_lines.append('*Contact pair, interaction=frictionless_contact')
    file_lines.append('\teval_tooth_neg.exposed_surface, Pulsator_jaw.y_min_surface')
    file_lines.append('*Contact pair, interaction=frictionless_contact')
    file_lines.append('\ttooth2_pos.exposed_surface, Pulsator_jaw.x_min_surface')

    for tooth in teeth:
        file_lines.append('*Boundary')
        file_lines.append('\t' + tooth.instance_name + '_pos.z0_nodes, 3, 3')
        file_lines.append('*Boundary')
        file_lines.append('\t' + tooth.instance_name + '_neg.z0_nodes, 3, 3')

    file_lines.append('*Boundary')
    file_lines.append('\tpulsator_jaw.z0_nodes, 3, 3')

    file_lines.append('*Boundary')
    file_lines.append('\tpinned, 1, 1')

    file_lines.append('*Boundary')
    file_lines.append('\t' + teeth[0].instance_name + '_pos.x1_nodes, 2, 2')
    file_lines.append('*Boundary')
    file_lines.append('\t' + teeth[-1].instance_name + '_neg.x1_nodes, 2, 2')

    file_lines.append('*Boundary')
    file_lines.append('\tjaw_ref_node, 3, 6')

    initiate_contact_lines = write_load_step('Initiate_contact', initial_inc=1e-4, force=1.)
    initiate_contact_lines.insert(3, '\t*Controls, reset')
    initiate_contact_lines.insert(4, '\t*Controls, parameters=line search')
    initiate_contact_lines.insert(5, '\t\t5, , , , ')
    initiate_contact_lines.insert(6, '\t*Contact Controls, Stabilize')
    initiate_contact_lines.insert(7, '\t*Contact Interference, shrink')
    initiate_contact_lines.insert(8, '\t\teval_tooth_neg.exposed_surface, Pulsator_jaw.y_min_surface')
    initiate_contact_lines.insert(9, '\t*Contact Interference, shrink')
    initiate_contact_lines.insert(10, '\t\ttooth2_pos.exposed_surface, Pulsator_jaw.x_min_surface')

    file_lines += initiate_contact_lines

    file_lines += write_load_step('loading', force=max_load*1e3, step_time=1.3, output_time_interval=0.02)

    with open(simulation_directory / 'pulsator_simulation.inp', 'w') as input_file:
        for line in file_lines:
            input_file.write(line + '\n')
        input_file.write("**EOF")


if __name__ == '__main__':
    main()
