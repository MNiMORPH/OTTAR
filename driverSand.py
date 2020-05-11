import riverwidth as rw
import importlib
importlib.reload(rw)

dlw = rw.DetachmentLimitedWidth(h_banks=1., S=1E-4, tau_crit=2., k_d=4E-6,
                                lambda_r=.1, b0=50.)

import numpy as np
t = np.arange(0,24*60.*60.*900,24*60.*60.)
Q = 200.*np.ones(len(t))

dlw.initialize(t,Q)
dlw.run()
dlw.finalize()
dlw.plot()
