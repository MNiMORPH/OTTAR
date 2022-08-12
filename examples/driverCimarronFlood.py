#! /usr/bin/python3

# Evolution of a river with cohesive banks

import ottar
from matplotlib import pyplot as plt

# For the sake of modifying the library: can reload and re-run
# Not necessary as part of an example
import importlib
importlib.reload(ottar)

# Create discharge time series
# In this case, a series of days that all have 200 m^3/s
import numpy as np
#"""
t = np.arange(0, 24*60.*60.*40, 24*60.*60./20.)
#Q = np.hstack((120000/3.28**3.*np.ones(5), 23000/3.28**3*np.ones(len(t)-255), 100*np.ones(250)))
Q = np.hstack((120000/3.28**3.*np.ones(5), 23000/3.28**3*np.ones(len(t)-755), 100*np.ones(750)))

# Instantiate river-width class
# Old
#rw = ottar.RiverWidth(h_banks=2., S=5E-4, tau_crit=2., k_d=4E-6,
#                                     b0=50., k_n=1.)
# New
rw = ottar.RiverWidth(h_banks=2., S=5E-4, tau_crit=.1, k_d=4E-5,
                                f_stickiness=.1, k_n_noncohesive=0.,
                                b0=50/3.28, D=1E-31)
#"""

"""
t = np.arange(0, 24*60.*60.*2000, 24*60.*60.)
Q = np.hstack((500.*np.ones(len(t))))#-250), 50*np.ones(250)))
#Q = np.hstack((500.*np.ones(250), 50*np.ones(len(t)-250)))
"""
# Instantiate river-width class
"""
# Old
rw = ottar.RiverWidth(h_banks=2., S=5E-4, tau_crit=2., k_d=4E-7,
                                     b0=50., k_n=.1)
# New
rw = ottar.RiverWidth(h_banks=2., S=5E-4, tau_crit=2, k_d=4E-7,
                                f_stickiness=.1, k_n_noncohesive=0.,
                                b0=50., D=1E-6)
"""


#t = t[:20]
#Q = Q[:20]

# Run full code: flow, widening, and plotting at the end
# rw.df is a Pandas DataFrame that holds the model outputs
rw.initialize_flow_calculations(0.03, 100, 1.5)
rw.initialize_timeseries(t,Q)
rw.run()
rw.finalize()
rw.plotQb()

