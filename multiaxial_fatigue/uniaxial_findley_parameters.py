import numpy as np

from scipy.optimize import fmin


def findley_stress(sa, sm, k):
    return 0.5*(k*(sa + sm) + np.sqrt(sa**2 + k**2*(sm + sa)**2))


def residual(parameters, datasets):
    residuals = [findley_stress(dataset[0], dataset[1], parameters[0]) - parameters[1] for dataset in datasets]
    return np.sum(np.array(residuals)**2)


def main():
    smooth_data = [(433, 90), (797, -297)]
    notched_data = [(440, 76), (825, -363)]
    k, sf = fmin(residual, [1., 600], args=(notched_data, ))
    print(findley_stress(smooth_data[0][0], smooth_data[0][1], k))
    print(findley_stress(smooth_data[1][0], smooth_data[1][1], k))
    print(sf)


if __name__ == '__main__':
    main()
