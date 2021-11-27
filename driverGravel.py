import riverwidth
from matplotlib import pyplot as plt
import importlib
importlib.reload(riverwidth)

rw = riverwidth.WidthNoncohesiveBanks(h_banks=4., S=1E-2, D=2E-2, Q0=5., k_n=1E-3)

import numpy as np
t = np.arange(0,10000*120,10000)
Q = 10.*np.ones(len(t))

rw.initialize(t,Q)
rw.run()
rw.finalize()
rw.plot()

