from hardess_convertion_functions import HV2HRC, HRC2HV

import numpy as np

hrc1 = HV2HRC(835)
hrc2 = HV2HRC(900)

au1 = 0.17
au2 = 0.1

print(np.linalg.solve(np.array([[au1, 1-au1], [au2, 1-au2]]), [hrc1, hrc2]))
print(HRC2HV(np.linalg.solve(np.array([[au1, 1-au1], [au2, 1-au2]]), [hrc1, hrc2])))
print(HRC2HV(35*0.1))
HRC = 63
dHv = (223*(100-HRC) + (223*HRC + 14500))/(100 - HRC)**2*35

print(0.05*dHv)
