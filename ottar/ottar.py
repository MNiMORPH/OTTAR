#! /usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import warnings
import yaml

class RiverWidth(object):
    """
    Transient adjustments to river-channel width
    """

    def __init__(self, yamlparams=None,
                       h_banks=None, S=None, b0=None, tau_crit=None,
                       k_d=None, k_E=None,
                       f_stickiness=0., k_n_noncohesive=0.,
                       Parker_epsilon=0.2, intermittency=1.,
                       D=None, rho_s=2650.,
                       tau_star_crit_sed=0.0495):

        # Make yamlparams (or lack thereof) a class variable
        self.yamlparams = yamlparams
    
        if yamlparams is not None:
            # streamflow
            self.streamflow_filename = str( yamlparams['streamflow']
                    ['filename'] )
            self.streamflow_datetime_column_name = str( yamlparams['streamflow']
                    ['datetime_column_name'] )
            self.streamflow_discharge_column_name = str( yamlparams['streamflow']
                    ['discharge_column_name'] )
            
            # widthdata
            # (optional, for comparison of model & data)
            self.widthdata_filename = str( yamlparams['widthdata']
                    ['filename'] )
            self.widthdata_datetime_column_name = str( yamlparams['widthdata']
                    ['datetime_column_name'] )
            self.widthdata_discharge_column_name = str( yamlparams['widthdata']
                    ['discharge_column_name'] )
            
            # morphology
            self.S = float( yamlparams['morphology']['slope'] )
            self.h_banks = float( yamlparams['morphology']['bank_height'] )
            # Input variable as initial state in channel-width list
            self.b = [ float( yamlparams['morphology']['initial_width'] ) ]
            self.bi = float( self.b[-1] )
            # back to morphology
            self.Parker_epsilon = float( 
                                  yamlparams['morphology']['Parker_epsilon'] )
            
            # sediment
            self.D = float( yamlparams['sediment']['grain_diameter'] )
            self.rho_s = float( yamlparams['sediment']['sediment_density'] )
            
            # widening
            self.tau_star_crit_sed = float( yamlparams['widening'] 
                  ['critical_shields_stress_of_noncohesive_sediment'] )
            self.tau_crit = float( yamlparams['widening'] 
                  ['critical_detachment_shear_stress_for_cohesive_sediement'] )
            self.k_E = float( yamlparams['widening'] 
                  ['noncohesive_entrainment_coefficient__k_E'] )
            self.k_d = float( yamlparams['widening'] 
                  ['coehsive_detachment_coefficient__k_d'] )
            
            # narrowing
            self.k_n_noncohesive = float( yamlparams['narrowing'] 
                    ['trapping_and_holding_efficiency__k_n_noncohesive'] )
            self.f_stickiness = float( yamlparams['narrowing'] 
                    ['f_stickiness'] )
            
            # doublemanning-flow
            self.manning2__n_ch = float( yamlparams['doublemanning-flow']
                    ['n_ch'] )
            self.manning2__k_fp = float( yamlparams['doublemanning-flow']
                    ['k_fp'] )
            self.manning2__P_fp = float( yamlparams['doublemanning-flow']
                    ['P_fp'] )
            self.manning2__stage_offset = float( 
                    yamlparams['doublemanning-flow']['stage_at_Q=0'] )

            
            # output
            

            # VARIABLES THAT ARE NOT SET
            self.intermittency = 1. # I might just remove this; save confusion
            
            
        else:
            if (h_banks is None) or (S is None) or (b0 is None):
                print("\nAt muminum, h_banks, S, and b0 must be set.\n")
                sys.exit(2)

            # Starting values (None) -- D in this too, but set in inputs
            self.channel_n = None
            self.tau_crit_sed = None
            #self.u_star_crit_sed = None
            
            # Input variables
            self.h_banks = h_banks
            self.S = S
            self.tau_crit = tau_crit # Critical stress to detach particles from bank
            self.k_d = k_d # Cohesive substrate: detachment-rate coefficient
            self.k_E = k_E # Noncohesive substrate: entrainment-rate coefficient
            self.f_stickiness = f_stickiness # Cohesive bank "stickiness": fraction
                                             # of sediment transferred to bank that
                                             # stays there & leads to narrowing 
            self.k_n_noncohesive = k_n_noncohesive # narrowing coefficient for
                                                   # noncohesive seds:
                                                   # lateral velocity efficiency
                                                   # (how much eddy diffusion
                                                   # affects the bed?) + sticking
            self.Parker_epsilon = Parker_epsilon
            self.intermittency = intermittency
            self.D = D # Grain diameter [m]

            # Input variable as initial state in list
            self.b = [b0]
            self.bi = self.b[-1]

            # For sediment (used in bed load calculations)
            self.tau_star_crit_sed = tau_star_crit_sed # Default: Wong & Parker 06

            self.rho_s = rho_s # Quartz density by default

        #############
        # CONSTANTS #
        #############
        self.g = 9.807
        self.rho = 1000.
        self.porosity = 0.35
        # For sediment (used in bed load calculations)
        self.SECONDS_IN_DAY = 86400

        #####################
        # DERIVED CONSTANTS #
        #####################
        if self.D is not None:
            self.tau_crit_sed = self.tau_star_crit_sed * \
                                  ( (self.rho_s - self.rho) * self.g * self.D)
            self.u_star_crit_sed = (self.tau_crit_sed / self.rho)**.5

        # Derived: Equilibrium width set by cohesion or grain size
        if self.tau_crit_sed is None:
            self.equilibrium_width_set_by_cohesion = True
        elif self.tau_crit is not None:
            if self.tau_crit >= self.tau_crit_sed:
                self.equilibrium_width_set_by_cohesion = True
        else:
            self.equilibrium_width_set_by_cohesion = False

        # Initialize list for all calculated flow depths
        self.h_series = [np.nan]
        
        # Initialize lists for the rates of widening and narrowing
        self.db_widening_cohesive_control = [np.nan] # F if D, T if cohesion
        self.db_widening_series = [np.nan]
        self.db_narrowing_series = [np.nan]
        
        # And a series for the shear stress on the banks
        # (Equals 1/(1+self.Parker_epsilon) * bed shear stress)
        self.tau_bank_series = [np.nan]
        
    @classmethod
    def from_yaml(cls, filepath):
        # "cls" is the class itself!
        # Thanks to Eric Hutton for sharing this magic
        with open(filepath) as fp:
            config = yaml.safe_load(fp)
        return cls(config)

    def dynamic_time_step(self, max_fract_to_equilib=0.1):
        # Currently part of a big, messy "update" step
        pass

    def initialize(self):
        """
        Using the YAML params file:
        * set up flow calculations
        * import time-series data
        """
        self.initialize_flow_calculations()
        self.initialize_timeseries()

    def initialize_flow_calculations(self, channel_n=None, fp_k=None, fp_P=None,
                                            stage_offset=None, use_Rh=None ):
        """
        Hard-code for double Manning
        """
        #if self.yamlparams:
        self.channel_n = channel_n
        self.hclass = FlowDepthDoubleManning(use_Rh)
        self.hclass.initialize( channel_n, fp_k, fp_P, stage_offset,
                                self.h_banks, self.b[-1], self.S)


    def initialize_timeseries(self, t=None, Q=None):
        self.t = list(t)
        self.Q = list(Q)

    def get_equilibriumWidth(self, Q_eq):
        """
        Steady-state width under erosion only as t-->infinity
        Only in-channel flow: Using provided h as bankfull
        No form drag assumed
        """
        
        # Cohesion sets ultimate width (Dunne & Jerolmack, 2018, & following)
        if self.equilibrium_width_set_by_cohesion:
            # Manning's n required
            if self.channel_n is None:
                warnings.warn("Manning's n must be set to compute the "+
                               "equilibrium width\n"+
                               "of a bank-cohesion-set channel.")
            # Ensure that this is still consistent with Nilay's updates to the
            # sand-bed equation set once she finishes
            else:
                self.k_b__eq = ( self.rho * self.g /
                                 ((1+self.Parker_epsilon) * 
                                 self.tau_crit) )**(5/3.) \
                                 * self.channel_n
                b_eq = self.k_b__eq * Q_eq * self.S**(7/6.)

        # Grain weight sets ultimate width (Parker, 1978, & following)
        else:
            self.k_b__eq = 0.17 / (self.g**.5 * 
                                    ((self.rho_s - self.rho)/self.rho)**(5/3.)
                                    * (1+self.Parker_epsilon)**(5/3.)
                                    * self.tau_star_crit_sed**(5/3.) )
            b_eq = self.k_b__eq * Q_eq * self.S**(7/6.) / self.D**1.5
        
        return b_eq

    def widen(self):
        """
        Compute river widening and append this to the list of river-widening
        rates.
        
        Use the slower of the two between cohesive and noncohesive as a
        rate-limiting step. It will be up to the user to input, for example,
        a very low cohesion (or 0) for something like a braided river with
        minimal organic/cohesive matter.
        """

        # Inefficient but complete: compute both and compare.
        # That is, unless D is not set, in which case, a fully cohesive system
        # is assumed
        
        # If not specified to evolve via bed load:
        # No grain size and/or no erosional coefficient
        if (self.D is None) or (self.k_E is None):
            self.db_widening = self.widen_cohesive()
            self.db_widening_cohesive_control.append( True )
        # Else if not specified to evolve via suspended load
        elif self.tau_crit is None:
            self.db_widening = self.widen_noncohesive()
            self.db_widening_cohesive_control.append( False )
        # Otherwise, use rate-limiting one of the both
        else:
            db_noncohesive = self.widen_noncohesive()
            db_cohesive = self.widen_cohesive()
            # Allow natural and smooth transition between which factor among
            # these controls the widening rate.
            # With a gravel-bed channel with gravel banks packed with mud
            # in between with a decent amount of cohesion, this should result
            # in the mud dominating the widening rate up until the very end,
            # when this is taken over by the difficulty of moving the clasts.
            # Record which one won in the list
            if db_noncohesive > db_cohesive:
                self.db_widening_cohesive_control.append( True )
                self.db_widening = db_cohesive
            else:
                self.db_widening_cohesive_control.append( False )
                self.db_widening = db_noncohesive

        # Record this into the list
        self.db_widening_series.append( self.db_widening )
        
        

    def widen_cohesive(self):
        """
        Widening set by shear stress on the banks and a linear (in stress) rate
        of particle detachment
        """
        if self.tau_bank > self.tau_crit:
            # 2* because erosion & deposition are symmetrical across banks
            if self.h < self.h_banks:
                db_widening = 2*self.k_d*self.h/self.h_banks \
                               * ( self.tau_bank - self.tau_crit ) \
                               * self.dt * self.intermittency
            else:
                db_widening = 2*self.k_d \
                               * ( self.tau_bank - self.tau_crit ) \
                               * self.dt * self.intermittency
        else:
            db_widening = 0.
        return db_widening

    def widen_noncohesive(self):
        """
        Widening set by sediment entrainment that is proportional to
        (tau* - tau*_c).
        All entrained sediment tumbles to the channel bed, where tau* is
        higher and the sediment is transported away.
        """

        if self.tau_star_bank > self.tau_star_crit_sed:
            # We assume the bank shape of Parker (1978), which is plausible
            # and allows the entire bank region to experience
            # the exact same shear stress regardless of flow depth.
            # We unexactly but not so unrealistically assume that this
            # shape is similar regardless of the flow depth.
            
            # 2* because erosion & deposition are symmetrical across banks

            # Incorporate symmetry, porosity, and (if a full hydrograph is
            # not being used) intermittency
            # Note: Porosity here but not in cohesive case (via convention
            # of empirical shear--detachment rule used there)
            
            # Overall general form (h<=h_banks):
            # 3.6 K h/h_banks (tau*-tau*_c) D/t_flight
            # 3.6 --> pi/4 * 4.6: coefficient for bed-sediment concentration
            # K --> coefficient; how many hops down bank? unconsidered terms?
            # h/h_banks --> length of erosion compared to total length to erode
            # D --> Erosional distance scale from one grain being removed
            # t_flight --> Length of time in transport; scales relationship
            #              between concentration and detachment rate
            #              (longer flight time --> lower detachment rate
            #              required to attain concentration)

            # WHAT TO DO IF CONC > 1?
            
            t_flight = 10.6 * (self.D / 
                                ( (self.rho_s - self.rho)/self.rho)
                                * self.g )**.5
            if self.h < self.h_banks:
                db_widening = 2 * 3.6 * self.k_E * self.h/self.h_banks \
                           * ( self.tau_star_bank - self.tau_star_crit_sed ) \
                           * self.D \
                           / t_flight \
                           * self.dt * self.intermittency \
                           / (1 - self.porosity)
            else:
                db_widening = 2 * 3.6 * self.k_E \
                           * ( self.tau_star_bank - self.tau_star_crit_sed ) \
                           * self.D \
                           / t_flight \
                           * self.dt * self.intermittency \
                           / (1 - self.porosity)
        # Otherwise, no erosion
        else:
            db_widening = 0.

        return db_widening
        
    def narrow_suspended_load(self):
        # Narrowing -- flipped sign convention of conc diff (so no $-$ needed)
        qsy = self.f_stickiness * self.K_Ey * \
                self.sed_conc_diff__suspended_load()
        return 2*qsy*self.dt / ( (1-self.porosity) * self.h_banks )

    def narrow_bed_load(self):
        """
        Lateral sediment velocity is proportional to downstream sediment
        velocity. This equation gives the bank-ward depth-integrated flux
        and then rescales it to bank lateral position change
        """
        if self.D is None:
            return 0
        if self.tau_star_bed < self.tau_star_crit_sed:
            return 0
        # Narrowing
        # Overall: bankward - channel-ward
        # qsby = u_{s,b}/4 * 3.6 (tau*_b-tau^*_c)/2 * (2/3)D
        # Channel to bank
        usx_ch = 4.4 * (self.u_star_bed - self.u_star_crit_sed) \
                    + 0.11 * ((self.rho_s - self.rho)/self.rho)**.5 \
                    * self.g**.5 * self.D**.5
        f_Am_ch = np.min( 
                    ( 3.6 * (self.tau_star_bed - self.tau_star_crit_sed),
                    1.) )
        # f_Am_ch/2. is because half of the sediment goes in each direction
        qsy_ch = usx_ch/4. * f_Am_ch/2. * 2/3.*self.D
        # Bank to channel
        if self.u_star_bank < self.u_star_crit_sed:
            qsy_bank = 0
        else:
            usx_bank = 4.4 * (self.u_star_bank - self.u_star_crit_sed) \
                        + 0.11 * ((self.rho_s - self.rho)/self.rho)**.5 \
                        * self.g**.5 * self.D**.5
            f_Am_bank = np.min( 
                            ( 3.6 * (self.tau_star_bank - self.tau_star_crit_sed),
                            1.) )
            qsy_bank = usx_bank/4. * f_Am_bank/2. * 2/3.*self.D
        # Together, net bankward transport
        qsy = qsy_ch - qsy_bank
        # Return effect on banks.
        # Multiplied by two: Both sides
        # Amplified by porosity in banks
        # Reduced by bank height (conversion to lateral rate of motion)
        return self.k_n_noncohesive \
                  *2*qsy*self.dt / ( (1-self.porosity) * self.h_banks )

    def compute__u_star__tau__K_Ey(self):
        self.u_star_bank = (self.tau_bank / self.rho)**.5
        if self.D is not None:
            self.tau_star_bank = self.tau_bank / ( (self.rho_s - self.rho) *
                                                    self.g * self.D )
        self.tau_bed = self.tau_bank * (1 + self.Parker_epsilon)
        self.u_star_bed = (self.tau_bed / self.rho)**.5
        if self.D is not None:
            self.tau_star_bed = self.tau_bed / \
                            ( (self.rho_s - self.rho) * self.g * self.D )
        self.K_Ey = 0.13 * self.h * self.u_star_bed

    def narrow(self):
        """
        Narrow based turbulent diffusion of sediment towards the banks
        (easy to visualize for suspended load; more of a discrete Brownian
        process for bed load)
        
        Result is the sum of supended-load and bed-load processes
        """

        if self.bi < 0:
            print("ERROR: Channel width < 0!")
            self.bi = 0
            raise ValueError('Negative flow depth given: Nonphysical.')            
        if self.h == 0:
            # No water in channel (!)
            print("Note: No water in channel.")
            self.db_narrowing = 0
        else:
            # Shear velocities and bed (center) shear stress
            # A bit of redundancy lies within.
            # 
            # K_Ey is the lateral eddy diffusivity [m^2/s].
            # Constant 0.13 is from Parker (1978, sand-bed)
            # Constant 0.16 (not used) was found in the work of Deng et al.
            # (2003) "Predicting Transverse Turbulent Diffusivity in Straight
            # Alluvial Rivers"
            # This is probably assuming that h < (b/2) or something
            # 
            # K_Ey Should also scale with lateral velocity variability that
            # moves bedload from side to side
            # This could be good to revisit experimentally
            self.K_Ey = 0.13 * self.h * self.u_star_bed

            # Noncohesive (bed load) + cohesive (suspended load)
            self.db_narrowing = self.narrow_bed_load() + \
                                self.narrow_suspended_load()
            
        # Record narrowing rate
        self.db_narrowing_series.append( self.db_narrowing )

        """
        if (self.D is not None) and (self.D > 0.002):
            _sus_load_narrowing = False
            sed_conc_diff_prop = self.sed_conc_diff__bed_load()
        else:
            _sus_load_narrowing = True
            sed_conc_diff_prop = self.sed_conc_diff__suspended_load()
        
        if sed_conc_diff_prop == 0:
            # Record narrowing rate
            self.db_narrowing = 0
            self.db_narrowing_series.append( self.db_narrowing )
            return # exit the function, returning None

        # Otherwise, continue.
            
        # Avoid div/0 (& math if not needed b/c no gradient)
        # As in, if there is somehow no water in the channel *and* sediment
        # is moving? Probably shouldn't need this. Raise warning.
        if sed_conc_diff_prop > 0:
            sed_conc_grad = sed_conc_diff_prop / min( self.h, self.h_banks )
        else:
            # Probably <0 (definitely?) because this would be sed motion
            # <=0 water
            warnings.warn("Water level <= 0 and almost certainly <0.")

        # K_Ey is the lateral eddy diffusivity [m^2/s].
        # Constant 0.13 is from Parker (1978, sand-bed)
        # Constant 0.16 (not used) was found in the work of Deng et al. (2003)
        # "Predicting Transverse Turbulent Diffusivity in Straight Alluvial
        # Rivers"
        # Should also scale with lateral velocity variability that will move
        # bedload from side to side -- though perhaps something more explict
        # would be good for bedload.
        # This is probably assuming that h < (b/2) or something
        K_Ey = 0.13 * self.h * self.u_star_bed
        
        # f_stickiness [unitless]: efficiency scaling term for lateral sediment
        #                          trapping and sticking

        # Lateral sediment discharge per unit channel length and width
        # [amount of sediment moving laterally / time]
        # (will declare it to be volumetric and define f_stickiness accordingly)
        if _sus_load_narrowing:
            self.qsy = self.f_stickiness * K_Ey * self.h * sed_conc_grad
        else:
            self.qsy = self.k_n_noncohesive * K_Ey * self.h * sed_conc_grad
        print(self.qsy)

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
        # Record narrowing rate
        self.db_narrowing_series.append( self.db_narrowing )
        return # Unnecessary but to make sure that the fucntion always returns
               # None (same type and value)
        """


    def sed_conc_diff__suspended_load(self):
        # Sediment concentrations / exit early if no change needed
        sed_conc_center_prop = (self.u_star_bed)**3.5
        if sed_conc_center_prop == 0:
            return 0 # exit the function, returning no gradient
        # Flows in excess of bank height are approximately unaffected by
        # differential drag: gradient tops out at that for h = h_banks
        if self.h < self.h_banks:
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
            
            # sed_conc_diff_prop
            return (sed_conc_center_prop - sed_conc_edge_prop)

        else:
            # u_star_banks up to h_banks, and then u_star_bed beyond
            # Could simplify
            return ( 1 - (1/(1+self.Parker_epsilon))**(7/4.) ) \
                    * (self.g * self.h_banks * self.S)**(3.5/2.) \

        # (later)
        # Concentration gradient calcualted over min(h, h_beta):
        # This sets distance from banks over which side-wall drag is important

    def update(self, dt, Qi, max_fract_to_equilib=0.1):
        """
        Euler forward wtih dynamic inner-loop time stepping
        Only widening; no narrowing
        """
        dt_outer = dt
        bi_outer = self.b[-1]
        self.hclass.set_b( bi_outer )
        h = self.hclass.compute_depth( Qi )
        # Use depth, not hydraulic radius, because the Parker_eposilon
        # factor is intended to convert the channel-centerline stress 
        # into a near-bank stress.
        self.tau_bank = self.rho * self.g * h * self.S / (1 + self.Parker_epsilon)
        # Record tau_bank in its series
        # Record narrowing rate
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
        # Update variables needed for computations
        self.compute__u_star__tau__K_Ey()
        # Compute widening
        self.widen()
        #self.b.append(self.bi + self.db_widening)
        self.narrow()
        self.h_series.append(h) # h is based on previous b but associated with
                                # the discharge that created current b
        self.b.append(self.bi + self.db_widening - self.db_narrowing)
        #print(self.hclass.compute_depth( 500. ))
        self.tau_bank_series.append( self.tau_bank )

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
        
        ##########################################
        # FINAL DATA PROCESSING AND ORGANIZATION #
        ##########################################
        
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

        # Generate additional numpy arrays
        self.db_widening_series = np.array(self.db_widening_series)
        self.db_narrowing_series = np.array(self.db_narrowing_series)
        
        # And create rates per day
        
        # Query if it is a timestamp
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            dt_ts = np.diff(self.t)
            dt_days = []
            for _dt in dt_ts:
                dt_days.append( _dt.total_seconds() / self.SECONDS_IN_DAY )
            dt_days = np.array(dt_days)
        # Assume numeric otherwise; maybe bad idea?
        else:
            dt_days = np.diff(self.t) / self.SECONDS_IN_DAY
        # nan at start -- don't have an initial elapsed time amount
        dt_days = np.hstack(( [np.nan], dt_days ))
        self.db_dt__day__widening_series = self.db_widening_series / dt_days
        self.db_dt__day__narrowing_series = self.db_narrowing_series / dt_days
        
        if self.yamlparams is not None:
            if self.yamlparams['output']['output_csv_filename'] is not None:
                self.write_csv( self.yamlparams['output']
                                               ['output_csv_filename'] )

        ############
        # PLOTTING #
        ############

        # If a YAML file is used, check whether plots should be generated
        # Empty lines default to False when recast from NoneType to Boolean
        if self.yamlparams is not None:
            if bool( self.yamlparams['plotting']
                  ['width'] ):
                self.plotb()
            if bool( self.yamlparams['plotting']
                  ['discharge+width'] ):
                self.plotQb()
            if bool( self.yamlparams['plotting']
                  ['widening+narrowing+stress'] ):
                self.plotWideningNarrowingStress()
            if bool( self.yamlparams['plotting']
                  ['width+widening+narrowing+stress'] ):
                self.plotWidthWideningNarrowingStress()
            if bool( self.yamlparams['plotting']
                  ['discharge+width+widening+narrowing+stress'] ):
                self.plotDischargeWidthWideningNarrowingStress()
            if bool( self.yamlparams['plotting']
                  ['discharge+width+widening+narrowing+grain-stress-ratio'] ):
                self.plotDischargeWidthWideningNarrowingGrainstressratio()
        
            # Then check whether they should be saved and/or shown
            # Currently, only one figure may be saved per run
            if self.yamlparams['plotting']['saveas'] is not None:
                plt.savefig( self.yamlparams['plotting']['saveas'] )
                
            # Then see if we should show the figure(s)
            if bool( self.yamlparams['plotting']['show'] ):
                plt.show()
                

    def plotb(self):
        """
        Plot channel width over time
        """
        # Manage typing of time stamp for plotting
        # Required with newer Pandas
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            _t = self.t
        else:
            _t = list(np.array(self.t)/86400.)
        #b_eq = self.get_equilibriumWidth(self.Qi)
        plt.figure()
        #plt.hlines(b_eq, self.t[0], self.t[-1]/86400.,
        #           '.5', label='Equilibrium width', linewidth=2)
        plt.plot(_t, self.b, 'k-', label='Transient width',
                 linewidth=2)
        plt.xlabel('Time [days of flood]')
        plt.ylabel('Channel width [m]')
        #plt.legend(loc='lower right')
        plt.tight_layout()
        # Show here by default if not running with a YAML configfile
        if self.yamlparams is None:
            plt.show()
        
    def plotQb(self, tdata=None, bdata=None):
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
        if tdata is not None and bdata is not None:
            ax1.plot(tdata, bdata, 'o', color='0.5')
        ax2.plot(_t, self.Q)
        ax2.set_ylabel('Discharge [m$^3$ s$^{-1}$]', fontsize=16)
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            ax2.set_xlabel('Date', fontsize=16)
        else:
            ax2.set_xlabel('Days since start', fontsize=16)
        plt.tight_layout()
        if self.yamlparams is None:
            plt.show()
        
    def plotWideningNarrowingStress(self, legend_loc=None):
        """
        Plot rates of widening and narrowing alongside bank stress
        """
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            _t = self.t
        else:
            _t = list(np.array(self.t)/86400.)
        plt.figure(figsize=(12,8))
        ax1 = plt.subplot(2,1,1)
        ax2 = plt.subplot(2,1,2)
        ax1.plot(_t, self.db_dt__day__widening_series, 'k-', linewidth=2, label='Widening')
        ax1.plot(_t, -self.db_dt__day__narrowing_series, '-', color='.5', linewidth=2, label='Narrowing')
        ax1.legend(loc=legend_loc, fontsize=12)
        ax1.set_ylabel('Channel width\nchange rate [m/day]', fontsize=16)
        ax2.plot(_t, self.tau_bank_series, 'k-', linewidth=2)
        ax2.set_ylabel('Bank shear stress [Pa]', fontsize=16)
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            ax2.set_xlabel('Date', fontsize=16)
        else:
            ax2.set_xlabel('Days since start', fontsize=16)
        plt.tight_layout()
        if self.yamlparams is None:
            plt.show()

        
    def plotWidthWideningNarrowingStress(self, legend_loc=None):
        """
        Plot width and rates of widening and narrowing alongside bank stress
        """
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            _t = self.t
        else:
            _t = list(np.array(self.t)/86400.)
        plt.figure(figsize=(8,8))
        ax1 = plt.subplot(3,1,1)
        ax2 = plt.subplot(3,1,2)
        ax3 = plt.subplot(3,1,3)
        ax1.plot(_t, self.b, 'k-', linewidth=2)
        ax1.set_ylabel('Channel width [m]', fontsize=12)
        ax2.plot(_t, self.db_dt__day__widening_series, 'k-', linewidth=2, label='Widening')
        ax2.plot(_t, -self.db_dt__day__narrowing_series, '-', color='.5', linewidth=2, label='Narrowing')
        ax2.legend(loc=legend_loc, fontsize=12)
        ax2.set_ylabel('Channel width\nchange rate [m/day]', fontsize=12)
        ax3.plot(_t, self.tau_bank_series, 'k-', linewidth=2)
        ax3.plot( [_t[0], _t[-1]] , [self.tau_crit, self.tau_crit], '--', 
                   color='.5', linewidth=1)
        ax3.set_ylabel('Bank shear stress [Pa]', fontsize=12)
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            ax3.set_xlabel('Date', fontsize=16)
        else:
            ax3.set_xlabel('Days since start', fontsize=16)
        plt.tight_layout()
        if self.yamlparams is None:
            plt.show()

        
    def plotDischargeWidthWideningNarrowingStress(self, legend_loc=None):
        """
        Plot discharge, width, and rates of widening and narrowing alongside
        bank stress
        """
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            _t = self.t
        else:
            _t = list(np.array(self.t)/86400.)
        plt.figure(figsize=(8,10))
        ax0 = plt.subplot(4,1,1)
        ax1 = plt.subplot(4,1,2)
        ax2 = plt.subplot(4,1,3)
        ax3 = plt.subplot(4,1,4)
        ax0.plot(_t, self.Q, 'k-', linewidth=2)
        ax0.set_ylabel('River discharge [m$^3$/s]', fontsize=12)
        ax1.plot(_t, self.b, 'k-', linewidth=2)
        ax1.set_ylabel('Channel width [m]', fontsize=12)
        ax2.plot(_t, self.db_dt__day__widening_series, 'k-', linewidth=2, label='Widening')
        ax2.plot(_t, -self.db_dt__day__narrowing_series, '-', color='.5', linewidth=2, label='Narrowing')
        ax2.legend(loc=legend_loc, fontsize=12)
        ax2.set_ylabel('Channel width\nchange rate [m/day]', fontsize=12)
        ax3.plot(_t, self.tau_bank_series, 'k-', linewidth=2)
        ax3.plot( [_t[0], _t[-1]] , [self.tau_crit, self.tau_crit], '--', 
                   color='.5', linewidth=1)
        ax3.set_ylabel('Bank shear stress [Pa]', fontsize=12)
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            ax3.set_xlabel('Date', fontsize=16)
        else:
            ax3.set_xlabel('Days since start', fontsize=16)
        plt.tight_layout()
        if self.yamlparams is None:
            plt.show()

        
    def plotDischargeWidthWideningNarrowingGrainstressratio(self, legend_loc=None):
        """
        Plot discharge, width, and rates of widening and narrowing alongside
        bank stress divided by critical stress to move particles
        """
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            _t = self.t
        else:
            _t = list(np.array(self.t)/86400.)
        plt.figure(figsize=(8,10))
        ax0 = plt.subplot(4,1,1)
        ax1 = plt.subplot(4,1,2)
        ax2 = plt.subplot(4,1,3)
        ax3 = plt.subplot(4,1,4)
        ax0.plot(_t, self.Q, 'k-', linewidth=2)
        ax0.set_ylabel('River discharge [m$^3$/s]', fontsize=12)
        ax1.plot(_t, self.b, 'k-', linewidth=2)
        ax1.set_ylabel('Channel width [m]', fontsize=12)
        ax2.plot(_t, self.db_dt__day__widening_series, 'k-', linewidth=2, label='Widening')
        ax2.plot(_t, -self.db_dt__day__narrowing_series, '-', color='.5', linewidth=2, label='Narrowing')
        ax2.legend(loc=legend_loc, fontsize=12)
        ax2.set_ylabel('Channel width\nchange rate [m/day]', fontsize=12)
        ax3.plot(_t, np.array(self.tau_bank_series)/self.tau_crit_sed, 'k-', linewidth=2)
        ax3.plot( [_t[0], _t[-1]] , [1, 1], '--', 
                   color='.5', linewidth=1)
        ax3.set_ylabel('Ratio: Bank stress \n/ critical stress', fontsize=12)
        if type(self.t[0]) == pd._libs.tslibs.timestamps.Timestamp:
            ax3.set_xlabel('Date', fontsize=16)
        else:
            ax3.set_xlabel('Days since start', fontsize=16)
        plt.tight_layout()
        if self.yamlparams is None:
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
    def __init__(self, use_Rh=True):
        # Default to using hyraulic radius and not just depth
        self.use_Rh = use_Rh

    def set_n(self, _var):
        self.n = _var

    def set_k(self, _var):
        self.k = _var

    def set_P(self, _var):
        self.P = _var
        
    def set_stage_offset(self, _var):
        self.stage_offset = _var

    def set_h_bank(self, _var):
        self.h_bank = _var

    def set_b(self, _var):
        self.b = _var

    def set_S(self, _var):
        self.S = _var

    def set_Q(self, _var):
        self.Q = _var

    def flow_depth_from_Manning_discharge( self, stage ):
        # flow depth
        h = stage - self.stage_offset
        # Does the flow go overbank?
        ob = h > self.h_bank
        if self.use_Rh:
            _r = h*self.b / (2*h + self.b)
        else:
            _r = h
        return self.b/self.n * _r**(5/3.) * self.S**0.5 \
                  + ob*self.k*(h-self.h_bank)**(ob*self.P) - self.Q

    def compute_depth(self, Q=None):
        if Q is not None:
            self.Q = Q
        if Q == 0:
            return 0
        else:
            return fsolve( self.flow_depth_from_Manning_discharge, 1. )[0]

    def initialize(self, n, k, P, stage_offset, h_bank, b, S):
        self.set_n(n)
        self.set_k(k)
        self.set_P(P)
        self.set_stage_offset(stage_offset)
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

