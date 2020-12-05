from hardess_convertion_functions import HV2HRC, HRC2HV

import numpy as np

hrc1 = HV2HRC(835)
hrc2 = HV2HRC(900)

au1 = 0.17
au2 = 0.1

print(np.linalg.solve(np.array([[au1, 1-au1], [au2, 1-au2]]), [hrc1, hrc2]))
print(HRC2HV(0.17*35 + 0.83*71))
print(HRC2HV(0.2*35 + 0.80*67))
print(HRC2HV(0.18*35 + 0.80*67 + 0.02*71))
print(HRC2HV(0.1*35 + 0.9*71))