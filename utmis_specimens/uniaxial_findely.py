import numpy as np
from scipy.optimize import fmin


def uniaxal_findley_stress(smax, smin, k):
    sa = (smax - smin)/2
    sm = (smax + smin)/2
    return 0.5*(k*(sa + sm) + (sa**2 + k**2*(sa + sm)**2)**0.5)


def residual(par, smax_data, smin_data):
    k = par[0]
    sf = par[1]
    res = 0
    for smax, smin in zip(smax_data, smin_data):
        res += (uniaxal_findley_stress(smax, smin, k) - sf)**2
    return res


if __name__ == '__main__':
    smooth_max = np.array([471.198, 505.])
    smooth_min = np.array([-1133.06, -372.176])
    notched_max = np.array([470.627,  495.152])
    notched_min = [-1192.57, -397.71]
    sa = (smooth_max - smooth_min)/2
    sm = (smooth_max + smooth_min)/2
    print((sa[1] - sa[0])/(sm[0] - sm[1]))

    sa = (notched_max - notched_min)/2
    sm = (notched_max + notched_min)/2
    print((sa[1] - sa[0])/(sm[0] - sm[1]))

    print(fmin(residual, [1, 500], args=(notched_max, notched_min)))
    k, sf = fmin(residual, [1, 500], args=(smooth_max, smooth_min))
    print(k, sf)

