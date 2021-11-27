import riverwidth as rw
import importlib
import pandas as pd
from matplotlib import pyplot as plt

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

t = data['Timestamp']
# First quarter of the record, repeated
#Q = list(data['Discharge (cfs)'] * cfs_to_m3s)[:10251] * 4
# Last quarter of the record, repeated
Q = list(data['Discharge (cfs)'] * cfs_to_m3s)[10251*3:] * 4

dlw.initialize(t, Q)
dlw.run()
dlw.finalize()

plt.figure()
#plt.semilogy(t/seconds_in_day, dlw.b, 'k-', label='Transient width',
#         linewidth=2)
#plt.title('Hydrological forcing from first quarter of 1906--2020 record')
plt.title('Hydrological forcing from fourth quarter of 1906--2020 record')
plt.plot(t, dlw.b, 'k-', label='Transient width',
         linewidth=2)
plt.plot(widthData['Year'], widthData['Width [m]'], 'o', color='0.5',
            label='Mankato-Henderson reach width (D. Libby)')
plt.ylabel('Channel width [m]')
plt.legend()
plt.tight_layout()
plt.show()

