import numpy as np


class InputFileReader:
    def __init__(self):
        self.nodal_data = None
        self.elements = {}
        self.set_data = {'nset': {}, 'elset': {}}

    def read_input_file(self, model_filename):
        nodes = []
        elements = {}
        with open(model_filename) as full_model_file:
            lines = full_model_file.readlines()
        key_word = None
        key_word_data = None
        for line in lines:
            if not line.startswith('**'):
                if line.startswith('*'):
                    key_word_line = line.split(',')
                    key_word = (key_word_line[0][1:]).lower().rstrip()   # Convert to lower case and remove newlines etc
                    key_word_data = [word.strip() for word in key_word_line[1:]]
                else:  # We have a data_row
                    data = line.split(',')
                    data = [item.rstrip() for item in data]
                    if key_word == 'node':
                        nodes.append(data)
                    elif key_word == 'element':
                        element_type = key_word_data[0][5:].rstrip()
                        if element_type not in elements:
                            elements[element_type] = []
                        elements[element_type].append(data)
                    elif key_word[-3:] == 'set':
                        set_name = key_word_data[0].split('=')[1].rstrip()
                        if set_name not in self.set_data[key_word]:
                            self.set_data[key_word][set_name] = []
                        self.set_data[key_word][set_name] += [int(label) for label in data if label]

        self.nodal_data = np.zeros((len(nodes), len(nodes[0])))

        for i, node in enumerate(nodes):
            self.nodal_data[i, :] = node

        for element_type, data in elements.items():
            self.elements[element_type] = np.array(data, dtype=int)

    def write_geom_include_file(self, filename, simulation_type='Mechanical'):
        file_lines = ['*NODE, NSET=ALL_NODES']
        for node in self.nodal_data:
            node_string = '\t' + str(int(node[0])) + ', '
            for n_data in node[1:]:
                node_string = node_string + str(n_data) + ', '
            file_lines.append(node_string[:-2])
        for element_type, element_data in self.elements.items():
            if element_type.endswith('R'):
                element_type = element_type[:-1]
            e_type = element_type

            if simulation_type != 'Mechanical':
                if element_type[1] == 'P':
                    e_type = 'DC2D' + element_type[-1]
                elif element_type[1] == 'A':
                    e_type = 'DCAX' + element_type[-1]
                else:
                    e_type = 'DC3D' + element_type[-1]
            file_lines.append('*ELEMENT, TYPE=' + e_type + ', ELSET=ALL_ELEMENTS', )
            for element in element_data:
                element_string = [str(e) for e in element]
                element_string = ', '.join(element_string)
                file_lines.append('\t' + element_string)

        with open(filename, 'w') as inc_file:
            for line in file_lines:
                inc_file.write(line + '\n')
            inc_file.write('**EOF')

    def write_sets_file(self, filename, skip_prefix='_', str_to_remove_from_setname='',
                        surfaces_from_element_sets=None):
        file_lines = []

        def write_set_rows(data_to_write):
            data_line = '\t'
            counter = 0
            for item in data_to_write:
                data_line += str(int(item)) + ', '
                counter += 1
                if counter == 16 or item == data_to_write[-1]:
                    file_lines.append(data_line[:-2])
                    counter = 0
                    data_line = '\t'

        for set_type, set_data in self.set_data.items():
            for key, data in set_data.items():
                key = key.upper()
                key = key.replace(str_to_remove_from_setname.upper(), '')
                if not key.startswith(skip_prefix) and (key.lower() not in ['all_elements', 'all_nodes']):
                    file_lines.append(('*' + set_type + ', ' + set_type + '=' + key).upper())
                    write_set_rows(data)

        if surfaces_from_element_sets:
            for surface_name, element_set_name in surfaces_from_element_sets:
                file_lines.append('*SURFACE, TYPE = ELEMENT, NAME=' + surface_name + ', TRIM=YES')
                file_lines.append('\t' + element_set_name)

        with open(filename, 'w') as set_file:
            for line in file_lines:
                set_file.write(line + '\n')
            set_file.write('**EOF')


if __name__ == '__main__':
    directory = '../fatigue_specimens/utmis_notched/'
    reader = InputFileReader()
    surfaces = [('EXPOSED_SURFACE', 'EXPOSED_ELEMENTS')]
    reader.read_input_file(directory + 'utmis_notched_geo.inc')
    reader.write_geom_include_file(directory + 'utmis_notched_geo_carbon.inc', simulation_type='Carbon')
    reader.write_geom_include_file(directory + 'utmis_notched_geo_thermal.inc', simulation_type='Thermal')
    reader.write_geom_include_file(directory + 'utmis_notched_geo_mechanical.inc', simulation_type='Mechanical')
    reader.write_sets_file(directory + 'utmis_notched_geo_sets.inc', str_to_remove_from_setname='Specimen_',
                           surfaces_from_element_sets=surfaces)
