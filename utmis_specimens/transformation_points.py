import numpy as np

from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib

from transformation_fatigue.materials.materials import SS2506

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


def transformation_start_load(data_set, a1, a2, a3, Mss):
    Ms = SS2506.Ms_1 + SS2506.Ms_2*data_set[:, 1]
    k = SS2506.k_1 + SS2506.k_2*data_set[:, 1]
    fm = 1 - np.exp(-k*(Ms + m_stress(data_set, a1, a2, a3) + Mss - 22.))
    r = fm - (1 - data_set[0, 2])
    return np.interp(0, r, data_set[:, 0])


def start_load_residual(parameters, data_sets):
    a1 = abs(parameters[0])
    a3 = abs(parameters[1])
    Mss = -abs(parameters[2])

    # Find a2 so that notched is at st
    st = 450

    def notched_residual(par, notched_data):
        fm0 = notched_data[0, 2]
        notched_data = notched_data[notched_data[:, 0] == [st], :]
        Ms = SS2506.Ms_1 + SS2506.Ms_2*notched_data[:, 1]
        k = SS2506.k_1 + SS2506.k_2*notched_data[:, 1]
        fm = np.exp(-k*(Ms + m_stress(notched_data, a1, par, a3) + Mss - 22.))
        return (fm[0] - fm0)**2*1000

    a2 = fmin(notched_residual, [1e-4], args=(data_sets[2],))[0]
    s = np.array([transformation_start_load(data_sets[0], a1, a2, a3, Mss),
                  transformation_start_load(data_sets[1], a1, a2, a3, Mss),
                  transformation_start_load(data_sets[2], a1, a2, a3, Mss)])
    smooth_start_load = (s[0] + s[1])/2 + max(s[0:2])
    r = 0
    r += 1000*np.count_nonzero(s == 0)
    return smooth_start_load + r


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
    b1, b2, b3, = 0.035, 1e-4, 0
    bss = M - m_stress(notched_center, b1, b2, b3)[-1]
    # par = fmin(start_load_residual, [0.02, 0, -56], args=([smooth_center, smooth_edge, notched_center],))
    # print(par)
    for data_set, su, c in zip([smooth_center, smooth_edge, notched_center], [424*2, 424*2, 237*2], ['b', 'g', 'r']):
        Mss, a1, a2, a3 = -128, 0.04, 4e-4, 0.e-07,
        fm = calculate_fm(data_set, a1, a2, a3, Mss)
        plt.figure(0)
        plt.plot(data_set[:, 0]/su, fm, c, lw=2)

        Mss, a1, a2, a3 = -105, 0.02, 2e-4, 2.e-07,
        fm = calculate_fm(data_set, a1, a2, a3, Mss)
        plt.figure(0)
        plt.plot(data_set[:, 0]/su, fm, c + ':', lw=2)

        plt.figure(1)
        plt.plot(data_set[:, 0]/su, -3*data_set[:, 3], '-' + c, lw=2)
        plt.plot(data_set[:, 0]/su, data_set[:, 4], '--' + c, lw=2)
        plt.plot(data_set[:, 0]/su, data_set[:, 5], ':' + c, lw=2)

    plt.figure(0)
    plt.xlim(0, 1)

    plt.show()


if __name__ == '__main__':
    main()
