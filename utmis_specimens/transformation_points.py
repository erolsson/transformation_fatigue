import numpy as np

from transformation_fatigue.materials.materials import SS2506

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def calculate_fm(data_set, a1, a2, a3, Mss):
    Ms = SS2506.Ms_1 + SS2506.Ms_2*data_set[:, 1]
    k = SS2506.k_1 + SS2506.k_2*data_set[:, 1]
    fm = 1 - np.exp(-k*(Ms + m_stress(data_set, a1, a2, a3) + Mss - 22.))
    fm[fm <= 1 - data_set[0, 2]] = 1 - data_set[0, 2]
    return fm


def m_stress(data_set, a1, a2, a3):
    return -3*a1*data_set[:, 3] + a2*data_set[:, 4]**2/3 + 2./27*a3*data_set[:, 5]**3


def main():
    s = np.linspace(0, 1000, 1000)
    uniaxial = np.zeros((1000, 6))
    uniaxial[:, 0] = s
    uniaxial[:, 1] = 0.008
    uniaxial[:, 2] = 0.2
    uniaxial[:, 3] = -s/3
    uniaxial[:, 4] = s
    uniaxial[:, 5] = s
    smooth_center = np.genfromtxt('utmis_smooth_center.csv', delimiter=',', skip_header=1)
    smooth_edge = np.genfromtxt('utmis_smooth_edge.csv', delimiter=',', skip_header=1)
    notched_center = np.genfromtxt('utmis_notched_center.csv', delimiter=',', skip_header=1)
    c = notched_center[0, 1]
    M = -np.log(notched_center[0, 2])/(SS2506.k_1 + SS2506.k_2*c) - SS2506.Ms_1 - SS2506.Ms_2*c + 22.
    b1, b2, b3, = 0.02, 1e-4, 0
    bss = M - m_stress(notched_center, b1, b2, b3)[-1]
    print(bss)
    for data_set, su, c in zip([smooth_center, smooth_edge, notched_center], [424*2, 424*2, 237*2], ['b', 'g', 'r']):
        Mss, a1, a2, a3 = -56.05272438731849, 0.02, 1e-4, 0
        fm = calculate_fm(data_set, a1, a2, a3, Mss)
        plt.figure(0)
        plt.plot(data_set[:, 0]/su, fm, c, lw=2)
        fm = calculate_fm(data_set, a1, a2, a3, -67.74183011161615)
        plt.figure(0)
        plt.plot(data_set[:, 0]/su, fm, ':' + c, lw=2)
        plt.figure(2)
        plt.plot(data_set[:, 0]/su, m_stress(data_set, a1, a2, a3) + Mss, c, lw=2)

        fm = calculate_fm(data_set, b1, b2, b3, bss)
        plt.figure(0)
        plt.plot(data_set[:, 0]/su, fm, '--' + c, lw=2)
        plt.figure(2)
        plt.plot(data_set[:, 0]/su, m_stress(data_set, b1, b2, b3) + bss, '--' + c, lw=2)

        plt.figure(1)
        plt.plot(data_set[:, 0]/su, -3*data_set[:, 3], c, lw=2)
        plt.plot(data_set[:, 0]/su, data_set[:, 4], '--' + c, lw=2)
        plt.plot(data_set[:, 0]/su, data_set[:, 5], ':' + c, lw=2)

    plt.figure(0)
    plt.xlim(0, 1)

    plt.figure(1)
    plt.xlim(0, 1)

    plt.figure(100)
    plt.plot(uniaxial[:, 0], calculate_fm(uniaxial, b1, b2, b3, bss))
    plt.show()


if __name__ == '__main__':
    main()
