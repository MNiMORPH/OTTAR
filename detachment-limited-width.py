#! /usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt

# Simple Euler forward

# Input variables
Q = 100.
b = [20.]
S = 1E-4
h_b = 1
lambda_r = .1 # m, guess, roughness. Should vary with h?
              # (depth-limited dunes?)
              # But is of limited importance -- raised to
              # VERY low power
epsilon = 0.2
tau_crit = 2 # Pa
kd = 4E-6 # m^2 s / kg
intermittency = 1

# Constants
g = 9.805
rho = 1000.

# Derived variables
a1 = rho * g**.7 * S**.7 * lambda_r**.1 / (8.1**.6 * (1+epsilon))

# Starting values
t = [0.]
dt = 24*60.*60.
nt = 900

# Equilibrium width?
beq = rho**(5/3.) * g**(7/6.) * S**(7/6.) * lambda_r**(1/6.) \
      / (8.1 * (1 + epsilon)**(5/3.)) * Q / tau_crit**(5/3.)

# Tau*
#tau_star_bed = h * S / ( ((rho_s - rho)/rho) * D)
#tau_star_bank = tau_star_bed / 1.2

# Compute through time
for i in range(nt):
    bi = b[-1]
    tau_bank = a1 * (Q/bi)**(3/5.)
    if tau_bank > tau_crit:
        bi += kd/h_b * ( tau_bank - tau_crit ) \
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
plt.show()
