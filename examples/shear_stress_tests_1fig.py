#! /usr/bin/python3

# Comparing stress and widening or narrowing rate

import ottar
from matplotlib import pyplot as plt

# Create discharge time series
# In this case, a series of days that all have 200 m^3/s
import numpy as np
#"""
t = np.arange(0, 24*60.*60.+1, 24*60.*60.)

Qi_range = np.logspace(-1,2,50)

# Widening -- cohesive only
# Strange scoop: Error in code?
tau_bank = []
h = []
db = []
for Qi in Qi_range:
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

# Bank height approximation -- not needed though
#Q_bf = np.where(h-rw.h_banks == np.min(np.abs(h - rw.h_banks)))[0][0]
tau_b_bf = rw.rho * rw.g * rw.h_banks * rw.S
tau_banks_bf = tau_b_bf / (1 + rw.Parker_epsilon)

fig = plt.figure(figsize=(6,8))

ax = plt.subplot(4,1,1)
ax.plot(tau_bank, db, 'k-')
yl = ax.get_ylim()
xl = ax.get_xlim()
ax.vlines(rw.tau_crit, yl[0], yl[1], colors='0.5', linestyles='dashed')
ax.text(rw.tau_crit + np.abs(np.diff(xl))*0.01,
        yl[1] - np.abs(np.diff(yl))*0.1,
        horizontalalignment='left', verticalalignment='top',
        s=r'$\tau_{c,c}$')
ax.vlines(tau_banks_bf, yl[0], yl[1], colors='0.5', linestyles='dashed')
ax.text(tau_banks_bf + np.abs(np.diff(xl))*0.01,
        yl[1] - np.abs(np.diff(yl))*0.1,
        horizontalalignment='left', verticalalignment='top',
        s=r'$\tau_{\beta,\mathrm{bf}}$')
ax.set_ylabel("$\dot{b}$: Widening Rate\n(Cohesive)\n[m/day]")

# Widening -- noncohesive only
tau_bank = []
db = []
for Qi in Qi_range:
    Q = np.hstack(Qi*np.ones(2))
    rw = ottar.RiverWidth(h_banks=1., S=1E-3, tau_crit=None, k_d=1E-6,
                                    f_stickiness=0., k_n_noncohesive=0.,
                                    b0=20., D=8E-3)
    rw.initialize_flow_calculations(0.03, 100, 1.5)
    rw.initialize_timeseries(t,Q)
    rw.run()
    rw.finalize()
    tau_bank.append(rw.tau_bank)
    db.append(np.diff(rw.b))
    
tau_star_bank = np.array(tau_bank) / ( (rw.rho_s - rw.rho) * rw.g * rw.D )
tau_star_ratio = tau_star_bank / rw.tau_star_crit_sed

tau_crit = rw.tau_star_crit_sed * ( (rw.rho_s - rw.rho) * rw.g * rw.D )

#plt.plot(tau_star_ratio, db, 'k-')
ax = plt.subplot(4,1,2)
ax.plot(tau_bank, db, 'k-')
yl = ax.get_ylim()
xl = ax.get_xlim()
ax.vlines(tau_crit, yl[0], yl[1], colors='0.5', linestyles='dashed')
ax.text(tau_crit + np.abs(np.diff(xl))*0.01,
        yl[1] - np.abs(np.diff(yl))*0.1,
        horizontalalignment='left', verticalalignment='top',
        s=r'$\tau_{n,c}$')
ax.set_ylabel("$\dot{b}$: Widening Rate\n(Noncohesive)\n[m/day]")

# Narrowing -- suspended load
tau_bank = []
db = []
for Qi in Qi_range:
    Q = np.hstack(Qi*np.ones(2))
    rw = ottar.RiverWidth(h_banks=1., S=1E-3, tau_crit=1., k_d=0.,
                                    f_stickiness=1., k_n_noncohesive=0.,
                                    b0=20., D=8E-3)
    rw.initialize_flow_calculations(0.03, 100, 1.5)
    rw.initialize_timeseries(t,Q)
    rw.run()
    rw.finalize()
    tau_bank.append(rw.tau_bank)
    db.append(np.diff(rw.b))
    
tau_ratio = np.array(tau_bank) / rw.tau_crit

ax = plt.subplot(4,1,3)
ax.plot(tau_bank, db, 'k-')
yl = ax.get_ylim()
xl = ax.get_xlim()
ax.set_ylabel("$\dot{b}$: Narrowing Rate\n(Suspended Load)\n[m/day]")

# Narrowing -- bed load
tau_bank = []
db = []
for Qi in Qi_range:
    Q = np.hstack(Qi*np.ones(2))
    rw = ottar.RiverWidth(h_banks=1., S=1E-3, tau_crit=100., k_d=0,
                                    f_stickiness=0., k_n_noncohesive=1.,
                                    b0=20., D=8E-3)
    rw.initialize_flow_calculations(0.03, 100, 1.5)
    rw.initialize_timeseries(t,Q)
    rw.run()
    rw.finalize()
    tau_bank.append(rw.tau_bank)
    db.append(np.diff(rw.b))
    
tau_star_bank = np.array(tau_bank) / ( (rw.rho_s - rw.rho) * rw.g * rw.D )
tau_star_ratio = tau_star_bank / rw.tau_star_crit_sed

tau_crit = rw.tau_star_crit_sed * ( (rw.rho_s - rw.rho) * rw.g * rw.D )

ax = plt.subplot(4,1,4)
ax.plot(tau_bank, db, 'k-')
yl = ax.get_ylim()
xl = ax.get_xlim()
ax.vlines(tau_crit, yl[0], yl[1], colors='0.5', linestyles='dashed')
ax.text(tau_crit + np.abs(np.diff(xl))*0.01,
        yl[0] + np.abs(np.diff(yl))*0.1,
        horizontalalignment='left', verticalalignment='bottom',
        s=r'$\tau_{n,c}$')
ax.set_xlabel(r"$\tau_\beta$: Bank Shear Stress [Pa]")
ax.set_ylabel("$\dot{b}$: Narrowing Rate\n(Bed Load)\n[m/day]")

plt.tight_layout()

plt.savefig('shear_stress__widening__narrowing.pdf')

plt.show()

