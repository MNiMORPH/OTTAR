import riverwidth as rw
import importlib
import pandas as pd
import datetime
from matplotlib import pyplot as plt
import numpy as np

importlib.reload(rw)
plt.ion()
cfs_to_m3s = 0.0283168466
seconds_in_day = float(24 * 60 * 60)

data = pd.read_csv('MN_River_Mankato_28042001_1903-06-01_2020-05-10.csv')
data['Timestamp'] = data['Timestamp'].astype('datetime64[ns]')
data = data.sort_values('Timestamp', ignore_index=True)

widthData = pd.read_csv('DevonLibby_Table3-3_Reach3.csv')
widthData['Year'] = pd.to_datetime(widthData['Year'], format='%Y')
widthData = widthData.sort_values('Year', ignore_index=True)

dlw = rw.DetachmentLimitedWidth(h_banks=1., S=1E-4, tau_crit=7., k_d=4E-6,
                                lambda_r=.1, b0=55.)

#t = list(data['Timestamp'])*15
#Q = list(data['Discharge (cfs)'] * cfs_to_m3s)*15
#t = np.arange(0, seconds_in_day*len(Q), seconds_in_day)

t = list(data['Timestamp'])
ts_end = datetime.datetime.timestamp(t[-1])
t_list_adendum = t[10251*3:]
ts_adendum = []
for _datetime in t_list_adendum:
    ts_adendum.append(datetime.datetime.timestamp(_datetime))
ts_adendum = np.array(ts_adendum)
ts_adendum += ts_end
ts_adendum += seconds_in_day # start 1 day after
_datetime_adendum = []
for _ts in ts_adendum:
    _datetime_adendum.append(datetime.datetime.fromtimestamp(_ts))
t += _datetime_adendum
# First quarter of the record, repeated
#Q = list(data['Discharge (cfs)'] * cfs_to_m3s)[:10251] * 4
# Last quarter of the record, repeated
Q = list(data['Discharge (cfs)'] * cfs_to_m3s) \
    + list(data['Discharge (cfs)'] * cfs_to_m3s)[10251*3:]

dlw.initialize_flow_calculations(0.030857097939076584, 183.21859785867875, 1.4791019236638066)
dlw.initialize_timeseries(t, Q)
dlw.run()
dlw.finalize()

tplot = np.array(t, dtype=np.datetime64)

plt.figure()
#plt.semilogy(t/seconds_in_day, dlw.b, 'k-', label='Transient width',
#         linewidth=2)
#plt.title('Hydrological forcing from first quarter of 1906--2020 record')
plt.title('Hydrological forcing from fourth quarter of 1906--2020 record')
plt.plot(tplot, dlw.b, 'k-', label='Transient width',
         linewidth=2)
plt.plot(widthData['Year'], widthData['Width [m]'], 'o', color='0.5',
            label='Mankato-Henderson reach width (D. Libby)')
plt.ylabel('Channel width [m]')
plt.legend()
plt.tight_layout()
plt.show()

