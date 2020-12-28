from collections import namedtuple

from transformation_fatigue.multiaxial_fatigue.multiaxial_fatigue_criteria import Findley


SteelData = namedtuple('Steel_data', ['HV'])


class _SS2506Material:
    def __init__(self):
        self.name = 'SS2506'

        self.composition = {'C': 0.23, 'Si': 0.2, 'Mn': 0.97, 'P': 0.01, 'S': 0.033, 'Cr':  0.6, 'Ni': 0.36, 'Mo': 0.7,
                            'Cu': 0.16, 'Al': 0.27}

        self.fatigue_life_knee = 1e5
        self.fatigue_exponent_low = 5.
        self.fatigue_exponent_high = 28.
        self.mean_stress_sensitivity_parameters = (0.017, 8.27e-4)

    def mean_stress_sensitivity(self, steel_properties, multiaxial_criterion):
        if multiaxial_criterion == Findley:
            return (self.mean_stress_sensitivity_parameters[0]
                    + self.mean_stress_sensitivity_parameters[1]*steel_properties.HV)

    @staticmethod
    def weibull_stress(steel_properties, multiaxial_criterion):
        if multiaxial_criterion == Findley:
            return 158.7 + 0.481538*steel_properties.HV

    @staticmethod
    def weibull_threshold_stress(steel_properties, multiaxial_criterion):
        return 0

    @staticmethod
    def weibull_exponent(steel_properties, multiaxial_criterion):
        if multiaxial_criterion == Findley:
            return 11.5719e6/steel_properties.HV**2

    @staticmethod
    def critical_effective_stress(steel_properties, multiaxial_criterion):
        if multiaxial_criterion == Findley:
            return 197.75 + 0.56833*steel_properties.HV


SS2506 = _SS2506Material()
