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

streamflowData = pd.read_csv('MN_Jordan_daily.csv')
streamflowData['Timestamp'] = streamflowData['Timestamp'].astype('datetime64[ns]')
streamflowData = streamflowData.sort_values('Timestamp')

rw = ottar.RiverWidth.from_yaml('config.yaml')

t = streamflowData['Timestamp']
Q = streamflowData['Discharge (cfs)'] * cfs_to_m3s

rw.initialize_flow_calculations(  )

rw.initialize_timeseries(t, Q)
rw.run()
rw.finalize()

