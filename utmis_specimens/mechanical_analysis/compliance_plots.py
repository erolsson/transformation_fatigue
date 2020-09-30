import os
import pickle

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    specimen_loads = {'smooth': {-1.: [737., 774., 820.], 0.: [425., 440.]},
                      'notched': {-1.: [427., 450.], 0.: [225., 240., 255.]}}

    for i, (specimen, load_cases) in enumerate(specimen_loads.items()):
        data = np.genfromtxt("compliance_utmis_" + specimen + ".csv", delimiter=",")
        plt.plot(data[:, 1], data[:, 2])
        elastic_compliance = (data[30, 2] - data[0, 2])/(data[30, 1] - data[0, 1])
        print(elastic_compliance)
        for R, amplitudes in load_cases.items():
            for amp in amplitudes:
                mean_load = (1 + R)/(1 - R)*amp
                max_load = mean_load + amp

                min_load = mean_load - amp
                max_strain = np.interp(max_load, data[:, 1], data[:, 2])
                mean_strain = np.interp(mean_load, data[:, 1], data[:, 2])
                plt.plot(max_load, max_strain, '*')
                P0 = max_strain/elastic_compliance - max_load
                print(P0)
                plt.plot([P0, max_load], [0, max_strain])
    plt.show()


if __name__ == '__main__':
    main()
