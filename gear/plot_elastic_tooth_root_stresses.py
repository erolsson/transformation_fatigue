import pathlib

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from abaqus_python.abaqus_interface import ABQInterface
from find_tooth_path_points import root_path

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def calculate_normal_stress(s, n):
    return sum(s[:3]*n**2) + 2*(s[3]*n[0]*n[1] + s[4]*n[0]*n[2] + s[5]*n[1]*n[2])


def main():
    abq = ABQInterface("abq2018")
    tooth_odb_directory = pathlib.Path.home() / "scania_gear_analysis" / "odb_files" / "mechanical_analysis"
    tooth_odb_filename = tooth_odb_directory / "elastic_stresses.odb"
    frames = abq.get_frames(tooth_odb_filename)
    model_directory = pathlib.Path.home() / "python_projects" / "python_fatigue" / "planetary_gear" / "input_files"
    inp_file = model_directory / 'quarter_tooth_tilt2.inp'
    root_path_points = root_path(inp_file)
    root_path_points[:, 0] *= -1
    v = np.array([
        root_path_points[-2, 0] - root_path_points[-1, 0],
        root_path_points[-2, 1] - root_path_points[-1, 1],
        0
    ])
    v /= np.linalg.norm(v)
    n = np.array([v[1], -v[0], 0])

    # Transform the vector n to cylinder coordinates
    r = np.sqrt(root_path_points[-1, 0]**2 + root_path_points[-1, 1]**2)
    q = np.arctan2(root_path_points[-1, 1], root_path_points[-1, 0])
    a = np.array([
        [np.cos(q), -np.sin(q), 0],
        [np.sin(q), np.cos(q), 0],
        [0, 0, 1]
    ])
    n = np.dot(a, n)
    sn = np.zeros(len(frames))
    force = 130*np.linspace(0, 1., len(frames))
    for frame in frames:
        s = abq.get_tensor_from_path(odb_file_name=tooth_odb_filename, path_points=root_path_points,
                                     field_id='S', frame_number=frame)
        sn[frame] = calculate_normal_stress(s[-1, :], n)
        print("The stress at force {F} kN is {s} MPa".format(F=force[frame], s=sn[frame]))
    plt.plot(force, sn)
    plt.xlabel('Force [kN]', fontsize=24)
    plt.ylabel('Tooth root stress [MPa]', fontsize=24)
    plt.tight_layout()
    plt.savefig('elastic_tooth_root_stresses.png')
    plt.show()


if __name__ == '__main__':
    main()
