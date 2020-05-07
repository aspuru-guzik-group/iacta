import numpy as np
import xtb_utils
from scipy.signal import find_peaks

xtb = xtb_utils.xtb_driver()

outputs = []
at1 = 0
at2 = 2
for i in range(100):
    inputf = "GabesCatalyst/init/opt%4.4i.xyz"%i
    outputs += [xtb.wf(inputf)]
    print(inputf, outputs[-1]['wbo'][at1,at2])
    
R = np.sqrt([np.sum((o['positions'][at1,:] - o['positions'][at2,:])**2)
             for o in outputs])
orders = np.array([o['wbo'][at1,at2] for o in outputs])

dO = abs(np.diff(orders)/np.diff(R))
indices, _ = find_peaks(abs(dO), prominence=0.5)
mtd = indices + 1

