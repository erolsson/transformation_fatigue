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
    temp = 22.
    carbon = np.linspace(0., 1, 1000)
    ea = -1.146e-2 + 4.22e-3*carbon + 2.377e-5*temp
    em1 = -4.740e-3 + 1.531e-02*carbon + -9.253e-03*carbon**2 + 1.210e-5*temp -4.361e-06*carbon*temp + 2.238e-09*temp**2
    a3 = -0.00175278
    a2 = 0.00939213
    a1 = -0.00389193
    em2 = a1 + a2*carbon + a3*carbon**2 + 1.342e-5*temp -3.92e-06*carbon*temp + 2.1857e-09*temp**2
    dv1 = (em1 - ea)*3
    dv2 = (em2 - ea)*3
    print(np.polyfit(carbon/100, dv2, 2))
    print(np.polyfit([0.0052, 0.01], [2.20242836e-02, 0.00586463], 1))
    plt.plot(carbon, dv1)
    plt.plot(carbon, dv2)
    plt.show()


if __name__ == '__main__':
    main()
