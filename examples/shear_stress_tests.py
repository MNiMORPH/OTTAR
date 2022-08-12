#! /usr/bin/python3

# Comparing stress and widening or narrowing rate

import ottar
from matplotlib import pyplot as plt

# Create discharge time series
# In this case, a series of days that all have 200 m^3/s
import numpy as np
#"""
t = np.arange(0, 24*60.*60.+1, 24*60.*60.)


# Widening -- cohesive only
# Strange scoop: Error in code?
tau_bank = []
h = []
db = []
for Qi in np.logspace(-1,3,50):
    Q = np.hstack(Qi*np.ones(2))
    rw = ottar.RiverWidth(h_banks=1., S=1E-3, tau_crit=2., k_d=1E-6,
                                    f_stickiness=0., k_n_noncohesive=0.,
                                    b0=20., D=None)
    rw.initialize_flow_calculations(0.03, 100, 1.5)
    rw.initialize_timeseries(t,Q)
    rw.run()
    rw.finalize()
    tau_bank.append(rw.tau_bank)
    h.append(rw.h)
    db.append(np.diff(rw.b))
    #print( rw.db_narrowing )
    print 
    
tau_ratio = np.array(tau_bank) / rw.tau_crit
h = np.array(h)

plt.figure()
plt.plot(tau_ratio, db, 'k-')
plt.plot(tau_ratio, h/rw.h_banks, 'b-')
plt.xlabel(r"$\tau_\beta / \tau_c$")
plt.ylabel("$\dot{b}$")
plt.title("Widening: Cohesive")
#plt.show()

# Widening -- noncohesive only
tau_bank = []
db = []
for Qi in np.logspace(-1,3,50):
    Q = np.hstack(Qi*np.ones(2))
    rw = ottar.RiverWidth(h_banks=1., S=1E-3, tau_crit=None, k_d=1E-6,
                                    f_stickiness=0., k_n_noncohesive=0.,
                                    b0=20., D=1E-2)
    rw.initialize_flow_calculations(0.03, 100, 1.5)
    rw.initialize_timeseries(t,Q)
    rw.run()
    rw.finalize()
    tau_bank.append(rw.tau_bank)
    db.append(np.diff(rw.b))
    
tau_star_bank = np.array(tau_bank) / ( (rw.rho_s - rw.rho) * rw.g * rw.D )
tau_star_ratio = tau_star_bank / rw.tau_star_crit_sed

plt.figure()
plt.plot(tau_star_ratio, db, 'k-')
plt.xlabel(r"$\tau^*_\beta / \tau^*_c$")
plt.ylabel("$\dot{b}$")
plt.title("Widening: Noncohesive")

# Narrowing -- suspended load
tau_bank = []
db = []
for Qi in np.logspace(-1,3,50):
    Q = np.hstack(Qi*np.ones(2))
    rw = ottar.RiverWidth(h_banks=1., S=1E-3, tau_crit=1., k_d=0.,
                                    f_stickiness=1., k_n_noncohesive=0.,
                                    b0=20., D=1E-4)
    rw.initialize_flow_calculations(0.03, 100, 1.5)
    rw.initialize_timeseries(t,Q)
    rw.run()
    rw.finalize()
    tau_bank.append(rw.tau_bank)
    db.append(np.diff(rw.b))
    
tau_ratio = np.array(tau_bank) / rw.tau_crit

plt.figure()
plt.plot(tau_bank, db, 'k-')
plt.xlabel(r"$\tau_\beta$")
plt.ylabel("$\dot{b}$")
plt.title("Narrowing: Suspended load")

# Narrowing -- bed load
tau_bank = []
db = []
for Qi in np.logspace(-1,3,50):
    Q = np.hstack(Qi*np.ones(2))
    rw = ottar.RiverWidth(h_banks=1., S=1E-3, tau_crit=100., k_d=0,
                                    f_stickiness=0., k_n_noncohesive=1.,
                                    b0=20., D=1E-2)
    rw.initialize_flow_calculations(0.03, 100, 1.5)
    rw.initialize_timeseries(t,Q)
    rw.run()
    rw.finalize()
    tau_bank.append(rw.tau_bank)
    db.append(np.diff(rw.b))
    
tau_star_bank = np.array(tau_bank) / ( (rw.rho_s - rw.rho) * rw.g * rw.D )
tau_star_ratio = tau_star_bank / rw.tau_star_crit_sed

plt.figure()
plt.plot(tau_star_ratio, db, 'k-')
plt.xlabel(r"$\tau^*_\beta / \tau^*_c$")
plt.ylabel("$\dot{b}$")
plt.title("Narrowing: Bed load")

plt.show()

