#! /usr/bin/python3

# Shear stress across a channel cross section
# Based on Eq. 25 & environs from Parker (1978, gravel-bed)
# Then I'm going to see find the magnitude of lateral stresses.

# Started 31 OCT 2021 by ADW in Potsdam @ Aaron Bufe's place!

import numpy as np
from matplotlib import pyplot as plt
plt.ion()

# Half widths
b_flat_half = 16 # [m]
b_bank_half = 4 # [m]

# Half-channel geometry
# Solve every cm
dx = 0.01
y = np.arange(0,20+1E-6,dx)
z_bed = 0*y # Assuming bed center is at 0
# 2 meters gained over 4 meters lateral distance
#z_bed[y>=b_flat_half] = 2 * (np.exp(y[y>=b_flat_half] - b_flat_half)-1) / (np.exp(4)-1) # too sharp
z_bed[y>=b_flat_half] = (y[y>=b_flat_half] - b_flat_half)**2 / 8

#plt.plot(y, z_bed)


# Looks all right. Next step, the equation.


rho = 1000.
g = 9.805

S = 0.001

rho_s = 2650.
D = 0.01 # 10 mm

# Channel depth at center
# Start wtih simple/bankfull
h_c = 2
h = h_c - z_bed

# Roughness = ks
# D84 20 mm?
k_s = 0.07

R_k = h_c/k_s # Garry's normalized depth thing

psi_approx = (1/12. * np.log(30*R_k) - 5/72.) * (1 + 1 / (2*np.log(30*R_k) - 17/3.))

B_margin_ratio = b_bank_half / (b_flat_half + b_bank_half)

epsilon = (2*h_c/(2*b_bank_half))**2

# Do we mean this theta, or the bed lateral angle theta?
theta = 1/psi_approx**0.5 * (1-B_margin_ratio) / (B_margin_ratio * epsilon**0.5)

eta = y / b_bank_half - (1 - B_margin_ratio)/B_margin_ratio



# Oof, need to redo the bank shape!
# Do so before any y derivatives

mu = 0.84 # internal friction on banks
beta = 0.85 # efficiency of lift

h2 = h_c / (1 - mu*beta) * ( np.cos( eta*np.arccos(mu*beta)) - mu*beta )
h[y >= b_flat_half] = h2[y >= b_flat_half]


# Approx curvature on boundary using central difference
d2h_dy2 = np.zeros(y.shape)
d2h_dy2[1:-1] = ( h[:-2] - 2*h[1:-1] + h[2:] ) / dx**2
# Will just calculate at the needed point in the future. For now:
d2h_dy2_junction = float( d2h_dy2[y == b_flat_half] )

d2f_deta2__yA = b_bank_half / h_c * d2h_dy2_junction
gamma = -0.5 * d2f_deta2__yA

# NOT NEEDED -- until gravity portion(below) is used
z_bed_ext = np.zeros(len(y)+1)
z_bed_ext[1:-1] = (z_bed[:-1] + z_bed[1:])/2.
z_bed_ext[0] = 0
z_bed_ext[-1] = (np.max(y)+dx - b_flat_half)**2 / 8
# Hard coded parabola here
theta_angle = np.arctan( np.diff(z_bed_ext)/dx )

d__h_c = 1 - epsilon * ( gamma * (1 + 2*psi_approx) / ( np.sinh(theta) + np.cosh(theta) ) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)
#d__h_c = 1 - epsilon * ( gamma * (1 - 2*psi_approx) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)
#d__h_c = 1 - epsilon * ( gamma * (1 - 2*psi_approx) / ( np.sinh(theta_angle) + np.cosh(theta_angle) ) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)


#/ ( np.sinh(theta) * np.cosh(theta)


# I lost a factor of two somewhere!
# Need to track down, but for now:
d__h_c = 1 - 2 * epsilon * ( gamma * (1 + 2*psi_approx) / ( np.sinh(theta) + np.cosh(theta) ) ) * np.cosh( 1/psi_approx**0.5 * y/h_c)
z_bed_effective = h_c - d__h_c * h_c

tau_b_x = d__h_c * rho*g*h_c*S
tau_b_y = -np.diff(tau_b_x)/dx * np.sqrt(2)/2.
tau_b_x_mid = (tau_b_x[:-1] + tau_b_x[1:])/2.

tau_b_angle_deg = np.arctan(tau_b_y / -tau_b_x_mid) * 180/np.pi

# For plotting
y_mid = (y[:-1] + y[1:])/2.
plt.plot(y_mid, tau_b_angle_deg)


# Gravity pulling rock down
tau_y_body = (rho_s - rho)*g*D*np.sin(theta_angle)
tau_y_body_mid = (tau_y_body[:-1] + tau_y_body[1:])/2.
# But not accounting for normal force, friction, repose, etc.

plt.plot(y_mid, tau_b_y - tau_y_body_mid)




tau_shields = tau_b_x / ((rho_s-rho)*g*D)












