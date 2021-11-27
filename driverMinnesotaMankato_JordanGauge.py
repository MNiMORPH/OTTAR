#! /usr/bin/python3

import riverwidth as rw
#import importlib
import pandas as pd
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters

register_matplotlib_converters()

#importlib.reload(rw)
#plt.ion()
cfs_to_m3s = 0.0283168466
seconds_in_day = float(24 * 60 * 60)

data = pd.read_csv('MN_River_Mankato_28042001_1903-06-01_2020-05-10.csv')
data['Timestamp'] = data['Timestamp'].astype('datetime64[ns]')
data = data.sort_values('Timestamp')

widthData = pd.read_csv('DevonLibby_Table3-3_Reach3.csv')
widthData['Year'] = pd.to_datetime(widthData['Year'], format='%Y')
widthData = widthData.sort_values('Year')

dlw = rw.WidthCohesiveBanks(h_banks=1., S=1E-4, tau_crit=6, k_d=2E-6,
                                b0=60.)

t = data['Timestamp']
Q = data['Discharge (cfs)'] * cfs_to_m3s

dlw.initialize_flow_calculations(0.030857097939076584, 183.21859785867875, 1.4791019236638066)
dlw.initialize_timeseries(t, Q)
dlw.run()
dlw.finalize()

plt.figure()
#plt.semilogy(t/seconds_in_day, dlw.b, 'k-', label='Transient width',
#         linewidth=2)
plt.plot(t, dlw.b, 'k-', label='Transient width',
         linewidth=2)
plt.plot(widthData['Year'], widthData['Width [m]'], 'o', color='0.5',
            label='Mankato-Henderson reach width (D. Libby)')
plt.ylabel('Channel width [m]')
plt.legend()
plt.tight_layout()
plt.show()

