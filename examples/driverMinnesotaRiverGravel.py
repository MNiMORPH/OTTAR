#! /usr/bin/python3

"""
This is meant as a schematic example, using data from the Minnesota River at
Jordan (MN) stream gauge. Parameters are not calibrated to data, and are not
intended to fit real river-width measurements.
"""

import ottar
import pandas as pd
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()

cfs_to_m3s = 0.0283168466
seconds_in_day = float(24 * 60 * 60)

streamflowData = pd.read_csv('MN_Jordan_daily.csv')
streamflowData['Timestamp'] = streamflowData['Timestamp'].astype('datetime64[ns]')
streamflowData = streamflowData.sort_values('Timestamp')

rw = ottar.RiverWidth(h_banks=6., S=1E-4, tau_crit=5, k_d=2E-7,
                                f_stickiness=5E-2, k_n_noncohesive=1E-2,
                                b0=65., D=.5E-3)

t = streamflowData['Timestamp']
Q = streamflowData['Discharge (cfs)'] * cfs_to_m3s

rw.initialize_flow_calculations(0.031, 180., 1.5)

rw.initialize_timeseries(t, Q)
rw.run()
rw.finalize()

rw.plotQb()

