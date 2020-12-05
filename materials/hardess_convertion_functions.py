from __future__ import print_function, division

import numpy as np


def HRC2HV(HRC):
    return (223.*HRC+14500)/(100-HRC)


def HV2HRC(HV):
    return (100.*HV - 14500)/(HV+223)


if __name__ == '__main__':
    HRC_data = np.arange(45., 65., 1.)
    for HRC in HRC_data:
        print('HRC:', HRC, '-> HV:', HRC2HV(HRC))
