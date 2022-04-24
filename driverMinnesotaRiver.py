#! /usr/bin/python3

import riverwidth as rw
import pandas as pd
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()

cfs_to_m3s = 0.0283168466
seconds_in_day = float(24 * 60 * 60)

streamflowData = pd.read_csv('MN_Jordan_daily.csv')
streamflowData['Timestamp'] = streamflowData['Timestamp'].astype('datetime64[ns]')
streamflowData = streamflowData.sort_values('Timestamp')

dlw = rw.WidthCohesiveBanks(h_banks=6., S=1E-4, tau_crit=5, k_d=3E-8,
                                b0=65.)

t = streamflowData['Timestamp']
Q = streamflowData['Discharge (cfs)'] * cfs_to_m3s

dlw.initialize_flow_calculations(0.031, 180., 1.5)

dlw.initialize_timeseries(t, Q)
dlw.run()
dlw.finalize()


dlw.plotQb()
