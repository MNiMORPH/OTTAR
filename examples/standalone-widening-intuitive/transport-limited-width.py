#! /usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt

# Simple Euler forward

# Input variables
Q = 10.
b = [20.]
S = 1E-2
D = 2E-2
h_b = 4
intermittency = 1

# Constants
phi = 3.97
g = 9.805
rho_s = 2700.
rho = 1000.
tau_star_crit = 0.0495

# Derived variables
a1 = 1. / h_b
a2 = S**0.7 / ( 2.9 * (rho_s - rho)/rho * g**0.3 * D**0.9 )
kh = D**.1 / (2.9 * g**.3 * S**.3)

# Starting values
t = [0.]
dt = 10000
nt = 120

# Equilibrium width?
beq = 0.17 / ( g**.5 * ((rho_s - rho)/rho)**(5/3.) * 1.2**(5/3.)
               * tau_star_crit**(5/3.) ) * Q * S**(7/6.) / D**1.5

# Depth?
h = kh * (Q/b[-1])**0.6

# Tau*
tau_star_bed = h * S / ( ((rho_s - rho)/rho) * D)
tau_star_bank = tau_star_bed / 1.2

# Compute through time
for i in range(nt):
    bi = b[-1]
    tau_star_bank = a2 * (Q/bi)**(3/5.) / 1.2
    if tau_star_bank > tau_star_crit:
        bi += a1 * ( tau_star_bank - tau_star_crit )**(3/2.) \
                 * dt * intermittency
    else:
        b = beq
        break
    b.append(bi)
    t.append(t[-1] + dt)

t = np.array(t)
b = np.array(b)

plt.figure()
plt.hlines(beq, t[0] / (24.*60.*60.), t[-1] / (24.*60.*60.),
           '.5', label='Equilibrium width', linewidth=2)
plt.plot(t / (24.*60.*60.), b, 'k-', label='Transient width',
         linewidth=2)
plt.xlabel('Flood duration [days]')
plt.ylabel('Channel width [m]')
plt.legend()
plt.tight_layout()
