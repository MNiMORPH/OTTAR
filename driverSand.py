import riverwidth
import importlib
importlib.reload(riverwidth)

rw = riverwidth.WidthCohesiveBanks(h_banks=1., S=1E-4, tau_crit=2., k_d=4E-6,
                                     b0=50.)

import numpy as np
t = np.arange(0,24*60.*60.*900,24*60.*60.)
Q = 200.*np.ones(len(t))

#rw.initialize(t,Q)
rw.initialize_timeseries(t,Q)
rw.run()
rw.finalize()
#rw.plot()
