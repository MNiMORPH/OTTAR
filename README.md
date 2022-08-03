[![DOI](https://zenodo.org/badge/261265317.svg)](https://zenodo.org/badge/latestdoi/261265317)

# :otter: OTTAR

Ode To Transient Ancho de los Rivers

Transiently evolving river-channel width as a function of streambank properties, sediment in transport, and the hydrograph.

## Purpose

This model is designed to compute the rates of river-channel widening and narrowing based on changing hydrological regimes. It is currently designed for rivers with cohesive banks, with a critical shear stress for particle detachment and an erosion-rate coefficient.

## Installation

From PyPI:
```sh
pip install ottar
```

Locally, inside a clone of this git repository (the `-e` permits you to make local updates to the code and have them incorporated into the way that OTTAR runs):
```sh
pip install -e .
```

## Structure

OTTAR contains:

* The `RiverWidth` class, which contains methods to evolve the width of an alluvial river.
* The `FlowDepthDoubleManning` class, which is used to estimate flow depth from discharge, even with an evolving river-channel geometry.

## Examples

There's a [folder for these](https://github.com/MNiMORPH/OTTAR/tree/master/examples)!

## Model inputs and outputs

### Inputs

#### Key input parameters (RiverWidth)

| **Variable** 	| **Description**                                                                                                                                                                                                                                                                                                   	| **Typical value(s)**        	|
|------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|-----------------------------	|
| `h_banks`    	    | **Stream-bank height**. This is the thickness of material that must be removed for the river to widen by one unit lateral distance.                                                                                                                                                                                            	| 1-5 m                       	|
| `S`          	    | **Channel downstream-directed slope**. This is used to compute shear stresses and (if necessary) flow depth from water discharge.                                                                                                                                                                                              	| 10<sup>-3</sup>             	|
| `b0`         	    | **Initial width**. Starting width of a channel.                                                                                                                                                                                                                                                                                	| 1&ndash;1000 m                   	|
| `tau_crit`   	    | **Critical shear stress required to start eroding muddy banks**. At this stress, the flow begins to be able to detach particles. When set up to perform an **inversion** using data on river widening and past flows, this is one of two key parameters to be estimated for rivers with detachment-limited banks. 	| 1&ndash;10 Pa                     	|
| `k_d`        	    | **Erosion-rate coefficient**. This determines the rate of erosion as a function of shear stress above critical. When set up to perform an **inversion** using data on river widening and past flows, this is the other of two key parameters to be estimated.                         	| ~10<sup>-7</sup> m / (Pa s) 	|
| `f_stickiness`    | **Fraction of suspended-load particles contacting the bank that "stick" to it**. This modulates the turbulence-driven lateral-transport term and its impact on channel-narrowing rate, and comprises the abillity of banks to trap sediment and to hold it.                                                                                                                                                                                                                                  	| ~0.01-1              	|
| `k_n_noncohesive` | **Narrowing coefficient (noncohesive sediment)**. Trapping and holding efficiency in regards to noncohesive sediment; this may be due to deep pits between grains and/or other bank-rougness properties.                                                                                                                                                                                                                                                                                     	| ~0.01-1              	|
| `Parker_epsilon`  | **Excess bed shear-stress factor**. $\tau_b = (1+\epsilon) \tau_\beta$, where $\tau_b$ is bed shear stress and $\tau_\beta$ is bank shear stress.                                                                                                                                                                           	| 0.2                            	|
| `intermittency`  | **Intermittency**. Fraction of the time that the discharge given is active. This is always equal to 1 for a full hydrograph, and is $\leq$ 1 when a characteristic "geomorphically effective" discharge is considered. It can be thought of as a time-dialation factor.                                                        	| 10<sup>-3</sup>-1                	|
| `D`          	    | **Sediment median grain size**. This is the median size of the material both in transport and in the banks, and is required for bedload and/or noncohesive-sediment-dominated systems. It may also be specified for rivers dominated by susepended load and bank cohesion, though will likely play a more minor role in these.	| 0.0002-1 m                       	|

#### Key input data sets and parameters (FlowDepthDoubleManning)

*This step is used to compute flow depths from a discharge time series, and may be skipped if you already posess a time series of flow depth*

* Discharge time series
* Manning's n (channel)
* Roughness / topogrpahy coefficient (floodplains)
* Depth / topography exponent (floodplains)

### Outputs

This program outputs a time series of channel width, `b(t)`. It organizes this within a Pandas DataFrame that can also be exported using the `write_csv()` function within the `RiverWidth` class.

Plots can also be made of just river width (`plotb()`) or of discharge and river width (`plotQb`).
