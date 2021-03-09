# river-channel-width
Transient and steady-state river-channel width calculations for both transport- and detachment-limited banks.

## Purpose

This model is designed to compute the rate of river-channel widening based on changing hydrological regimes. This works for either:
* Detachment-limited banks: Our focus is cohesive muds, but this should, in principle, also work for rock or biological mats. Think rivers with muddy banks.
* Transport-limited banks. In this case, there is minimal to no vegetation or mud holding the banks together, or they do not add any appreciable strength. In this case, the rate of channel widening is based purely on the rate at which material (gravel, sand) can be removed from the banks and transported away.

## Model inputs and outputs

### Inputs

Required inputs from field data:
* River bank height, measured from the channel bed to the top of the banks
* Time-series of **flow depth** or a way to calculate this internally

Inputs, either from field data or to be calibrated in an inversion scheme (see descriptions on table below):
* Detachment-limited case
  * `tau_crit`
  * `k_d`
* Transport-limited case
  * `D`

#### Key input parameters

| **Variable** 	| **Description**                                                                                                                                                                                                                                                                                                   	| **Typical value(s)**        	|
|--------------	|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|-----------------------------	|
| `h_banks`    	| **Stream-bank height**. This is the thickness of material that must be removed for the river to widen by one unit lateral distance.                                                                                                                                                                               	| 1-5 m                       	|
| `S`          	| **Channel downstream-directed slope**. This is used to compute shear stresses and (if necessary) flow depth from water discharge.                                                                                                                                                                                 	| 10<sup>-3<\sup>             	|
| `tau_crit`   	| **Critical shear stress required to start eroding muddy banks**. At this stress, the flow begins to be able to detach particles. When set up to perform an **inversion** using data on river widening and past flows, this is one of two key parameters to be estimated for rivers with detachment-limited banks. 	| 1--2 Pa                     	|
| `k_d`        	| **Erosion-rate coefficient**. This determines the rate of erosion as a function of shear stress above critical. When set up to perform an **inversion** using data on river widening and past flows, this is the other of two key parameters to be estimated for rivers with detachment-limited banks.            	| 10<sup>--6<\sup> m / (Pa s) 	|
| `lambda_r`   	| **Nikuradse roughness coefficient** ($$k_s$$). This is used if required to create a stage--discharge relationship.                                                                                                                                                                                                	| 0.05                        	|
| `D`          	| **Bed-material median grain size**. This is required for rivers with transport-limited banks. When set up to perform an **inversion** using data on river widening and past flows, this is (if unknown) the parameters to be estimated for rivers with transport-limited banks.                                   	| 0.002--0.2 m                	|
| `Q0`         	| **Initial discharge**. Initial discharge to which a channel's width is adjusted                                                                                                                                                                                                                                   	| 1--1000 m<sup>3<\sup>/s     	|
| `b0`         	| **Initial width**. Starting width of a channel                                                                                                                                                                                                                                                                    	| 1--1000 m                   	|

### Outputs

This program outputs a time series of channel width, `b(t)`.
