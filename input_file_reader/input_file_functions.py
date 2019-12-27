import re

import numpy as np


def read_nodes_and_elements(model_filename):
    nodes = []
    elements = []
    with open(model_filename) as full_model_file:
        lines = full_model_file.readlines()
        reading_nodes = False
        reading_elements = False

        for line in lines:
            if bool(re.search('node', line.split()[0], re.IGNORECASE)):
                reading_nodes = True
                reading_elements = False
            elif bool(re.search('element', line.split()[0], re.IGNORECASE)):
                reading_elements = True
                reading_nodes = False
            elif '*' in line.split()[0]:
                reading_nodes = False
                reading_elements = False
            elif reading_nodes:
                nodes.append(line.split(","))
            elif reading_elements:
                elements.append(line.split(","))
    nodal_data = np.zeros((len(nodes), 4))

    for i, node in enumerate(nodes):
        nodal_data[i, :] = node
    elements = np.array(elements, dtype=int)
    return nodal_data, elements


def get_elements_from_nodes(node_labels, all_elements):
    nodal_id_set = set(node_labels)
    # Find the elements corresponding to the model nodes
    element_data = []
    for element in all_elements:
        include = True
        for node in element[1:]:
            if int(node) not in nodal_id_set:
                include = False
        if include:
            element_data.append([int(e) for e in element])
    return np.array(element_data, dtype=int)


def write_set_rows(data_to_write, file_lines):
    data_line = ''
    counter = 0
    for item in data_to_write:
        data_line += str(int(item)) + ', '
        counter += 1
        if counter == 16 or item == data_to_write[-1]:
            file_lines.append(data_line[:-2])
            counter = 0
            data_line = ''


def write_sets(node_sets, element_sets):
    file_lines = ['** Include file for sets in a quarter model of a planetary gear for dante sim']

    for key, data in element_sets.iteritems():
        file_lines.append('*Elset, elset=' + key)
        write_set_rows(data, file_lines)

    for key, data in node_sets.iteritems():
        file_lines.append('*Nset, nset=' + key)
        write_set_rows(data, file_lines)
    return file_lines


def write_geom_include_file(nodal_data, element_data, filename, simulation_type='Mechanical'):
    element_type = 'DC3D8'
    if simulation_type == 'Mechanical':
        element_type = 'C3D8'

    file_lines = ['*NODE']
    for node in nodal_data:
        file_lines.append('\t' + str(int(node[0])) + ', ' + str(node[1]) + ', ' + str(node[2]) + ', ' + str(node[3]))
    file_lines.append('*ELEMENT, TYPE=' + element_type)
    for element in element_data:
        element_string = [str(e) for e in element]
        element_string = ', '.join(element_string)
        file_lines.append('\t' + element_string)

    with open(filename, 'w') as inc_file:
        for line in file_lines:
            inc_file.write(line + '\n')
        inc_file.write('**EOF')
