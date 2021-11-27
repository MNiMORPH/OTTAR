import riverwidth as rw
from matplotlib import pyplot as plt
import importlib
importlib.reload(rw)

tlw = rw.WidthNoncohesiveBanks(h_banks=4., S=1E-2, D=2E-2, Q0=5., k_n=1E-3)

import numpy as np
t = np.arange(0,10000*120,10000)
Q = 10.*np.ones(len(t))

tlw.initialize(t,Q)
tlw.run()
tlw.finalize()

plt.figure()
plt.plot(t/(60*60*24.), tlw.b)
plt.xlabel('Days')
plt.ylabel('Width [m]')
plt.tight_layout()
plt.show()
