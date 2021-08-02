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
    smooth_data = np.array([
        [66.412, 438.58799999999997],
        [-330.931, 802.1289999999999]
    ])

    notched_data = np.array([
        [48.721000000000004, 446.431],
        [-360.9715, 831.5985]
    ])

    ks = -(smooth_data[0, 1] - smooth_data[1, 1])/(smooth_data[0, 0] - smooth_data[1, 0])
    kn = -(notched_data[0, 1] - notched_data[1, 1])/(notched_data[0, 0] - notched_data[1, 0])
    su_s = smooth_data[0, 1] + ks*smooth_data[0, 0]
    print(ks, su_s)
    su_n = notched_data[0, 1] + kn*notched_data[0, 0]
    print(kn, su_n)
    print(750/kn, 750/ks)

    plt.plot(smooth_data[:, 0], smooth_data[:, 1], 'rx', ms=12, mew=3)
    plt.plot(notched_data[:, 0], notched_data[:, 1], 'bx', ms=12, mew=3)
    plt.plot(380.5, 380.5, 'kx', ms=12, mew=3)
    plt.ylim(0)
    plt.show()


if __name__ == '__main__':
    main()
