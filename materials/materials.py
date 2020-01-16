from __future__ import print_function
import numbers

import numpy as np


# noinspection PyPep8Naming
class ElasticPlasticTransformMaterial:
    def __init__(self, E, v, sy0M, sy0A, Q, b, Cm, gamma_m, a, Ms, name, Mss, fM, beta, alpha, n, sde,
                 g0, g1, g2, g_mean, g_std, fsb0=0.):
        # Elastic parameters
        self.E = float(E)
        self.v = float(v)

        self.G = self.E/2/(1+self.v)
        self.K = self.E/3/(1-2*self.v)

        # Initial yield stress of Martensite and Austenite
        self.sy0M = sy0M
        self.sy0A = sy0A
        self.fM = fM
        # Parameters for isostropic hardening
        self.Q = Q
        self.b = b

        # parameters for Kinematic hardening
        self.Cm = Cm
        self.gamma_m = gamma_m

        # Martensite start temperature
        self.Ms = Ms

        # parameters for phase transformations
        self.a1 = a[0]
        self.a2 = a[1]
        self.a3 = a[2]

        # Parameters for the plastic strain transformations
        self.beta = beta
        self.alpha = alpha
        self.n = n

        self.R1 = 0.02
        self.R2 = 0.012

        self.dV = 0.037

        self.k = 0.017

        self.sde = sde

        self.g0 = g0
        self.g1 = g1
        self.g2 = g2

        self.g_mean = g_mean
        self.g_std = g_std

        self.fsb0 = fsb0

        self.M_sigma = 22
        self.M_d = 383.585

        self.back_stresses = Cm.shape[0]
        self.name = name
        if isinstance(Mss, numbers.Real):
            self.Mss = Mss
        else:
            self.Mss = (-1./self.k*np.log(1 - Mss[0]) - self.Ms - Mss[2]*(self.a1 + self.a2) +
                        Mss[1])

    def umat_depvar(self):
        return 15 + (self.back_stresses + 1)*6

    def material_input_file_string(self):
        parameters = [self.E, self.v, self.sy0M[0], self.sy0M[1], self.sy0A[0],  self.sy0A[1], self.Q[0], self.Q[1],
                      self.b[0], self.b[1], len(self.gamma_m)]
        kinematic_hardening_params = []
        for C, g in zip(self.Cm, self.gamma_m):
            kinematic_hardening_params += [C[0], C[1], g[0], g[1]]
        parameters += kinematic_hardening_params + [self.sde, self.R1, self.R2, self.dV, self.Ms, self.Mss,
                                                    self.k, self.a1, self.a2, self.a3, self.beta, self.alpha,
                                                    self.n, self.g0, self.g1, self.g2, self.g_mean, self.g_std,
                                                    self.M_sigma, self.M_d]

        file_lines = ['*Material, name=' + self.name,
                      '\t*Depvar',
                      '\t\t' + str(self.umat_depvar()),
                      '\t\t1, PLASTIC_STRAIN',
                      '\t\t2, AUSTENITE',
                      '\t\t3, CARBON',
                      '\t\t4, FERRITE',
                      '\t\t5, HARDNESS',
                      '\t\t6, LBAINITE',
                      '\t\t7, PEARLITE',
                      '\t\t8, Q_MARTENSITE',
                      '\t\t9, T_MARTENSTE',
                      '\t\t10, U_BAINITE',
                      '\t\t11, S_MARTENSITE',
                      '\t\t12, E_MARTENSITE',
                      '\t\t13, FSB',
                      '\t\t14, FSB0',
                      '\t\t15, R']
        counter = 16
        for m in range(len(self.Cm)):
            for comp in ['11', '22', '33', '12', '13', '23']:
                file_lines.append('\t\t' + str(counter) + ', BACK_STRESS_' + str(m) + '_' + comp)
                counter += 1
        for comp in ['11', '22', '33', '12', '13', '23']:
            file_lines.append('\t\t' + str(counter) + ', TOTAL_BACK_STRESS_' + comp)
            counter += 1

        file_lines.append('\t*User Material, constants=' + str(len(parameters)) + ', unsymm')
        parameter_str = ''
        for i, par in enumerate(parameters):
            if i % 8 == 0 and i != 0:
                file_lines.append('\t\t' + parameter_str)
                parameter_str = ''
            elif i != 0:
                parameter_str += ', '
            parameter_str += str(par)
        file_lines.append('\t\t' + parameter_str)
        return file_lines

    def sy0(self, fm):
        return fm*self.sy0M + (1 - fm)*self.sy0A


SS2506_no_trans = ElasticPlasticTransformMaterial(E=205e3, v=0.27, sy0M=(-662.15481215, 2.28816129), sy0A=(-253., 0.874),
                                                  Q=(0., 0.), b=(100., 0.),
                                                  Cm=np.array([(5767.48146841, 37.40347908),
                                                               (456688.92231459, -377.9659805),
                                                               (-397639.4939082, 1180.67400981)]),
                                                  gamma_m=np.array([(119.64906081, -0.10120371),
                                                                    (2077.94869114, -2.30167388),
                                                                    (19024.56244157, -23.72193205)]),
                                                  a=0*np.array([1., 0., 0.]),
                                                  Ms=20, name='SS2506_no_trans', Mss=-48.194154926470915, fM=0.78,
                                                  beta=0., alpha=115.6731692918583, n=4., sde=0.00,
                                                  g0=3.077651873747456,
                                                  g1=68.83381914607745, g2=0, g_mean=0, g_std=29.5540022577844,
                                                  fsb0=0.22758101605717712)

SS2506 = ElasticPlasticTransformMaterial(E=205e3, v=0.27, sy0M=(-662.15481215, 2.28816129), sy0A=(-253., 0.874),
                                         Q=(0., 0.), b=(100., 0.),
                                         Cm=np.array([(5767.48146841, 37.40347908),
                                                      (456688.92231459, -377.9659805),
                                                      (-397639.4939082, 1180.67400981)]),
                                         gamma_m=np.array([(119.64906081, -0.10120371),
                                                           (2077.94869114, -2.30167388),
                                                           (19024.56244157, -23.72193205)]),
                                         a=0.0*np.array([1., 1., 0.])/2,
                                         Ms=166.69903883197608, name='SS2506', Mss=-143.44, fM=0.8,
                                         beta=841.893, alpha=129.5, n=4., sde=0.00,
                                         g0=-1.918/2, g1=5.18, g2=1.918/2., g_mean=0, g_std=1.,
                                         fsb0=0.12948)

SS2506.k = 0.0132
SS2506.R1 = 0.035577
SS2506.R2 = 0.0
SS2506.dV = 0.0348
"""
a1=0.0020086602368106056, 
g_std=29.5540022577844, 
beta=306.6016112399211, 
R1=0.018685671366484635, 
R2=0.005105251428246899, 
fsb0=0.22758101605717712, 
dV=0.022435285414492856, 
alpha=115.6731692918583, 
g0=115.6731692918583, 
g1=68.83381914607745,
"""
# Hazar et. al
# Temperature   beta
# 22            ??
# 75            4.92067433e+02

if __name__ == '__main__':
    for line in SS2506_no_trans.material_input_file_string():
        print(line)
