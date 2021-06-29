import pathlib

import numpy as np

from input_file_reader.input_file_reader import InputFileReader


def get_tooth_profile(inp_file):
    reader = InputFileReader()
    reader.read_input_file(inp_file)
    exposed_nodes = set(reader.set_data['nset']['exposed_nodes'])
    nodes = reader.nodal_data[reader.nodal_data[:, 3] == 0]
    tooth_profile = []
    for node in nodes:
        if node[0] in exposed_nodes:
            tooth_profile.append([node[1], node[2]])
    tooth_profile = np.array(tooth_profile)
    ymin = np.min(tooth_profile[:, 1])
    ymax = np.max(tooth_profile[:, 1])
    tooth_profile = tooth_profile[tooth_profile[:, 1] > ymin + 0.2*(ymax - ymin), :]

    return tooth_profile[np.argsort(tooth_profile[:, 0]), :]


def root_path(input_file, tangent_angle=30., path_points=100, length=3.):
    tooth_profile = get_tooth_profile(input_file)
    x = tooth_profile[:, 0]
    y = tooth_profile[:, 1]
    x = x[y < (np.max(y) + np.min(y))/2]
    y = y[y < (np.max(y) +
               np.min(y))/2]
    tangent = np.diff(y[0] - y)/np.diff(x)
    angle = 180/np.pi*np.arctan(tangent)
    x_mid = x[0:-1] + np.diff(x)/2
    idx = np.argsort(angle)
    angle = angle[idx]
    x_mid = x_mid[idx]
    x_tangent = np.interp(90-tangent_angle, angle, x_mid)
    y_tangent = np.interp(x_tangent, x, y)
    path_data = np.zeros((path_points, 3))
    r = np.linspace(0, length, 100)
    path_data[:, 0] = r*np.cos(tangent_angle*np.pi/180)
    path_data[:, 1] = r*np.sin(tangent_angle*np.pi/180)
    path_data[:, 0] += -path_data[-1, 0] + x_tangent
    path_data[:, 1] += -path_data[-1, 1] + y_tangent
    path_data[:, 2] = 1e-3
    return path_data


def main():
    import matplotlib.pyplot as plt
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    inp_file = model_directory / 'quarter_tooth_tilt2.inp'
    tooth_profile = get_tooth_profile(inp_file)
    rp = root_path(inp_file)

    plt.figure(0)
    plt.plot(tooth_profile[:, 0], tooth_profile[:, 1], '-*')
    plt.plot(rp[:, 0], rp[:, 1], '-')
    plt.axis('scaled')

    plt.show()


if __name__ == '__main__':
    main()
