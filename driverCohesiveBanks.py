#! /usr/bin/python3

# Evolution of a river with cohesive banks

import riverwidth
from matplotlib import pyplot as plt

# For the sake of modifying the library: can reload and re-run
# Not necessary as part of an example
import importlib
importlib.reload(riverwidth)

# Create discharge time series
# In this case, a series of days that all have 200 m^3/s
import numpy as np
t = np.arange(0,24*60.*60.*1500,24*60.*60.)
Q = 200.*np.ones(len(t))

# Instantiate river-width class
rw = riverwidth.WidthCohesiveBanks(h_banks=2., S=1E-4, tau_crit=2., k_d=4E-6,
                                     b0=50., k_n=0)

# Run full code: flow, widening, and plotting at the end
# rw.df is a Pandas DataFrame that holds the model outputs
rw.initialize_flow_calculations(0.03, 100, 1.5)
rw.initialize_timeseries(t,Q)
rw.run()
rw.finalize()
rw.plot()

