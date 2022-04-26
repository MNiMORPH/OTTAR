#! /usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

class RiverWidth(object):
    """
    Transient adjustments to river-channel width
    """

    def __init__(self, h_banks, S, tau_crit, k_d, b0, k_n=0.,
                 Parker_epsilon=0.2, intermittency=1.):

        # Input variables
        self.h_banks = h_banks
        self.S = S
        self.tau_crit = tau_crit # Critical stress to detach particles from bank
        self.k_d = k_d # Bank rate constant
        self.intermittency = intermittency
        self.Parker_epsilon = Parker_epsilon
        self.k_n = k_n # Narrowing coefficient

        # Input variable as initial state in list
        self.b = [b0]
        self.bi = self.b[-1]

        # Constants
        self.g = 9.805
        self.rho = 1000.
        self.porosity = 0.35

        # Initialize list for all calculated flow depths
        self.h_series = [np.nan]

    def dynamic_time_step(self, max_fract_to_equilib=0.1):
        # Currently part of a big, messy "update" step
        pass

    def initialize_flow_calculations(self, channel_n, fp_k, fp_P ):
        """
        Hard-code for double Manning
        """
        self.hclass = FlowDepthDoubleManning()
        self.hclass.initialize( channel_n, fp_k, fp_P,
                                self.h_banks, self.b[-1], self.S)

        # Derived constant -- Manning's n required
        self.k_b__eq = ( self.rho * self.g /
                         ((1+self.Parker_epsilon) * self.tau_crit) )**(5/3.) \
                       * channel_n


    def initialize_timeseries(self, t, Q):
        self.t = list(t)
        self.Q = list(Q)

    def get_equilibriumWidth(self, Q_eq):
        """
        Steady-state width under erosion only as t-->infinity
        Only in-channel flow: Using this as bankfull
        No form drag assumed
        """
        b_eq = self.k_b__eq * Q_eq * self.S**(7/6.)
        return b_eq

    def widen(self):
        """
        Widen a river channel based on the shear stress (compared to a
        threshold) at the bank
        """
        if self.tau_bank > self.tau_crit:
            # 2* because erosion & deposition are symmetrical across banks
            if self.h < self.h_banks:
                self.db_widening = 2*self.k_d*self.h/self.h_banks \
                                     * ( self.tau_bank - self.tau_crit ) \
                                     * self.dt * self.intermittency
            else:
                self.db_widening = 2*self.k_d \
                                     * ( self.tau_bank - self.tau_crit ) \
                                     * self.dt * self.intermittency
        else:
            self.db_widening = 0.

    def narrow(self):
        """
        Narrow based turbulent diffusion of sediment towards the banks
        (easy to visualize for suspended load, more of a discrete Brownian
        process for bed load)

        Here for suspended load
        """

        # Shear velocities and bed (center) shear stress
        # A bit of redundancy lies within
        self.u_star_bank = (self.tau_bank / self.rho)**.5
        self.tau_bed = self.tau_bank * (1 + self.Parker_epsilon)
        self.u_star_bed = (self.tau_bed / self.rho)**.5

        # Sediment concentrations / exit early if no change needed
        sed_conc_center_prop = (self.u_star_bed)**3.5
        if sed_conc_center_prop == 0:
            self.db_narrowing = 0
            return # exit the function, returning None
        sed_conc_edge_prop = (self.u_star_bank)**3.5

        # Had previously divided by "bi" to create the gradient, but this
        # was forgetting my own work by hand! So much of the channel is
        # unaffected by the walls (constant velocity and stress), such that
        # the near-wall zone of substantital velocity and stress gradients
        # has a separate and approximately constant width.
        # Now, the only width feedback will be the related to the powers
        # to which the stress gradient is raised.
        # Realism & increased stability! (Though narrow channels still can have
        # deeper flows & steeper gradients)
        sed_conc_grad_prop = (sed_conc_center_prop - sed_conc_edge_prop)

        # Concentration gradient calcualted over min(h, h_beta):
        # This sets distance from banks over which side-wall drag is important

        # Avoid div/0 (& math if not needed b/c no gradient)
        if sed_conc_grad_prop > 0:
            sed_conc_grad = sed_conc_grad_prop / min( self.h, self.h_banks )

        # K_Ey is the lateral eddy diffusivity [m^2/s].
        # Constant 0.13 is from Parker (1978, sand-bed)
        # Constant 0.16 (not used) was found in the work of Deng et al. (2003)
        # "Predicting Transverse Turbulent Diffusivity in Straight Alluvial
        # Rivers"
        # This is probably assuming that h < (b/2) or something
        K_Ey = 0.13 * self.h * self.u_star_bed
        
        # k_n [unitless]: efficiency scaling term for lateral sediment
        #                 transport, trapping, and deposition

        # Lateral sediment discharge per unit channel length and width
        # [amount of sediment moving laterally / time]
        # (will declare it to be volumetric and define k_n accordingly)
        self.qsy = self.k_n * K_Ey * self.h * sed_conc_grad

        # EVENTUALLY HOLD TWO OPTIONS HERE:
        # 1. Uniform narrowing across the full h_banks
        # 2. Narrowing up to the height of the water only (so happens faster),
        #    and tracking the width of a virtual inset channel until the next
        #    highest flow moves the sediment farther up.
        #    even allowing the sediment to move up instantaneously will require
        #    some amount of tracking an arbitrary number of inset rectangles
        #    and therefore some extra coding and bookkeeping
        #    so let's include this, but not just yet.

        # SIMPLE NARROWING (OVER FULL CHANNEL WALL HEIGHT)
        # 2* because erosion & deposition are symmetrical across banks
        # Divided by amount of lateral sediment motion needed to create a unit
        # of narrowing across the full h_banks
        self.db_narrowing = 2*self.qsy*self.dt \
                                / ( (1-self.porosity) * self.h_banks )
        return # Unnecessary but to make sure that the fucntion always returns
               # None (same type and value)


    def update(self, dt, Qi, max_fract_to_equilib=0.1):
        """
        Euler forward wtih dynamic inner-loop time stepping
        Only widening; no narrowing
        """
        dt_outer = dt
        bi_outer = self.b[-1]
        self.hclass.set_b( bi_outer )
        h = self.hclass.compute_depth( Qi )
        self.tau_bank = self.rho * self.g * h * self.S / (1 + self.Parker_epsilon)
        if self.tau_bank > self.tau_crit:
            self.bi = self.b[-1]
            dt_remaining = dt
            while dt_remaining != 0:
                if dt_remaining < 0:
                    raise RuntimeError('More time used than allowed. '
                                          +str(dt_remaining)
                                          +' seconds remaining')
                self.tau_bank = self.rho * self.g * h * self.S / (1 + self.Parker_epsilon)
                dbdt = 2*self.k_d*h/self.h_banks \
                           * ( self.tau_bank - self.tau_crit ) \
                           * self.intermittency
                """
                b_eq = np.inf # self.get_equilibriumWidth(self.Qi)
                dt_to_cutoff = max_fract_to_equilib * (b_eq - self.bi) / dbdt
                dt_inner = np.min((dt_to_cutoff, dt_remaining))
                self.bi += self.k_d/self.hself.b[-1]_banks \
                              * ( self.tau_bank - self.tau_crit ) \
                              * dt_inner * self.intermittency
                dt_remaining -= dt_inner
                #print(dt_remaining, self.bi, b_eq)
                """
                self.bi += dbdt * dt_outer
                dt_remaining = 0 # Perhaps return inner loop later
        self.h_series.append(h) # h is based on previous b but associated with
                                # the discharge that created current b
        self.b.append(self.bi) # add current b to list of b

    def update__simple_time_step(self, dt, Qi):
        """
        Simple Euler forward.
        Has widening and narrowing.
        """
        self.bi = self.b[-1]
        # Update for the Manning calc
        self.hclass.set_b( self.bi )
        # Might find a different way for this going forward
        self.dt = dt
        # Current discharge and shear stress
        # Is this updated for the rating-curve 2x Manning approach?
        h = self.hclass.compute_depth( Qi )
        self.tau_bank = self.rho * self.g * h * self.S \
                / (1 + self.Parker_epsilon)
        self.h = h # For the widening, at least for now
        # Compute widening
        self.widen()
        #self.b.append(self.bi + self.db_widening)
        self.narrow()
        self.h_series.append(h) # h is based on previous b but associated with
                                # the discharge that created current b
        self.b.append(self.bi + self.db_widening - self.db_narrowing)
        #print(self.hclass.compute_depth( 500. ))

    def run(self):
        for i in range(1, len(self.t)):
            # Not sure how inefficient this repeat check will be
            # Will find a better way in the future
            try:
                dt = (self.t[i] - self.t[i-1]).total_seconds()
            except:
                dt = (self.t[i] - self.t[i-1])
            #self.update(dt, self.Q[i])
            self.update__simple_time_step(dt, self.Q[i])

    def finalize(self):
        """
        Prepending NaN at the beginning of discharge, water_depth, and time
        lists. This is because width list is initalized with b0. Water depth is
        prepended at the beginning of the class.
        
        The output DataFrame will contain:
        * Timestamp [datetime]
        * Discharge [m^3/s]
        * Channel width [m]
        * Water depth [m]
        
        This is (consider updating to make width, depth consistent):
        Timestamp[now]
        Discharge[now]
        Depth[from discharge[now] and width[before]]
        Width[[new]]

        """
        
        # Generate numpy arrays of equal length
        self.t = np.array(self.t)
        self.b = np.array(self.b)
        self.Q = np.array(self.Q)
        self.h_series = np.array(self.h_series)
        
        # Create Pandas dataframe
        self.df = pd.DataFrame()
        self.df['Timestamp'] = self.t
        self.df['Discharge [m^3/s]'] = self.Q
        self.df['Channel width [m]'] = self.b
        self.df['Water depth [m]'] = self.h_series

    def plotb(self):
        """
        Plot channel width over time
        """
        #b_eq = self.get_equilibriumWidth(self.Qi)
        plt.figure()
        #plt.hlines(b_eq, self.t[0], self.t[-1]/86400.,
        #           '.5', label='Equilibrium width', linewidth=2)
        plt.plot(self.t/86400., self.b, 'k-', label='Transient width',
                 linewidth=2)
        plt.xlabel('Time [days of flood]')
        plt.ylabel('Channel width [m]')
        #plt.legend(loc='lower right')
        plt.tight_layout()
        plt.show()
        
    def plotQb(self):
        """
        Plot channel width and river discharge over time
        """
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            _t = self.t
        else:
            _t = list(np.array(self.t)/86400.)
        plt.figure(figsize=(12,8))
        ax1 = plt.subplot(2,1,1)
        ax2 = plt.subplot(2,1,2)
        ax1.plot(_t, self.b, 'k-', linewidth=2)
        ax1.set_ylabel('Channel width [m]', fontsize=16)
        ax2.plot(_t, self.Q)
        ax2.set_ylabel('Discharge [m$^3$ s$^{-1}$]', fontsize=16)
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            ax2.set_xlabel('Date', fontsize=16)
        else:
            ax2.set_xlabel('Days since start', fontsize=16)
        plt.tight_layout()
        plt.show()
        
    def write_csv(self, filename):
        """
        Write CSV of the DataFrame:
        * Date [datetime]
        * Discharge [m^3/s]
        * Channel width [m]
        * Water depth [m]
        """
        if filename[-4:] != '.csv':
            filename += '.csv'
        self.df.to_csv(filename)


###############################
## FLOW DEPTH FROM DISCHARGE ##
###############################

from scipy.optimize import fsolve

class FlowDepthDoubleManning( object ):
    """
    Use Manning's equation to obtain flow depth from river discharge,
    using a conversion from ManningFit.py outputs
    """
    def __init__(self):
        pass

    def set_n(self, _var):
        self.n = _var

    def set_k(self, _var):
        self.k = _var

    def set_P(self, _var):
        self.P = _var

    def set_h_bank(self, _var):
        self.h_bank = _var

    def set_b(self, _var):
        self.b = _var

    def set_S(self, _var):
        self.S = _var

    def set_Q(self, _var):
        self.Q = _var

    def flow_depth_from_Manning_discharge( self, h ):
        ob = h > self.h_bank
        return self.b/self.n * h**(5/3.) * self.S**0.5 \
                  + ob*self.k*(h-self.h_bank)**(ob*self.P) - self.Q

    def compute_depth(self, Q=None):
        if Q is not None:
            self.Q = Q
        if Q == 0:
            return 0
        else:
            return fsolve( self.flow_depth_from_Manning_discharge, 1. )[0]

    def initialize(self, n, k, P, h_bank, b, S):
        self.set_n(n)
        self.set_k(k)
        self.set_P(P)
        self.set_h_bank(h_bank)
        self.set_b(b)
        self.set_S(S)

    def update(self, Q=None):
        """
        Not exactly updating anything, but to follow standard CSDMS I(U)RF
        """
        self.h = self.compute_depth(Q)
        return self.h

    def run(self, Q=None):
        """
        Not exactly running anything, but to follow standard CSDMS I(U)RF
        Same as "update" step
        """
        return self.update(Q)

    def finalize(self):
        pass
