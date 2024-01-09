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

rw = ottar.RiverWidth.from_yaml('config.yaml')

"""
rw = ottar.RiverWidth(h_banks=5.8, S=1E-4, tau_crit=5, k_d=2E-6, k_E=0.4,
                                b0=65., f_stickiness=1E-2,
                                k_n_noncohesive=2E-5, D=0.00025
                                )
"""

t = streamflowData['Timestamp']
Q = streamflowData['Discharge (cfs)'] * cfs_to_m3s

#t = streamflowData['Timestamp'][:4000]
#Q = streamflowData['Discharge (cfs)'][:4000] * cfs_to_m3s

rw.initialize_flow_calculations( 0.034, 138., 1.62,
                                  stage_offset=0.47, use_Rh=True )


rw.initialize_timeseries(t, Q)
rw.run()
rw.finalize()

rw.plotQb()

