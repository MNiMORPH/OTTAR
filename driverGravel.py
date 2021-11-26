import riverwidth as rw
import importlib
importlib.reload(rw)

tlw = rw.WidthNoncohesiveBanks(h_banks=4., S=1E-2, D=2E-2, Q0=5.)

import numpy as np
t = np.arange(0,10000*120,10000)
Q = 10.*np.ones(len(t))

tlw.initialize(t,Q)
tlw.run()
tlw.finalize()

