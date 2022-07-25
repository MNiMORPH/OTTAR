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
|--------------	|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|-----------------------------	|
| `h_banks`    	| **Stream-bank height**. This is the thickness of material that must be removed for the river to widen by one unit lateral distance.                                                                                                                                                                               	| 1-5 m                       	|
| `S`          	| **Channel downstream-directed slope**. This is used to compute shear stresses and (if necessary) flow depth from water discharge.                                                                                                                                                                                 	| 10<sup>-3</sup>             	|
| `tau_crit`   	| **Critical shear stress required to start eroding muddy banks**. At this stress, the flow begins to be able to detach particles. When set up to perform an **inversion** using data on river widening and past flows, this is one of two key parameters to be estimated for rivers with detachment-limited banks. 	| 1&ndash;10 Pa                     	|
| `k_d`        	| **Erosion-rate coefficient**. This determines the rate of erosion as a function of shear stress above critical. When set up to perform an **inversion** using data on river widening and past flows, this is the other of two key parameters to be estimated.            	| ~10<sup>-7</sup> m / (Pa s) 	|
| `k_n`         	| **Narrowing coefficient**. This modulates the efficiency of channel narrowing via lateral sediment transport and deposition. It may relate to bar/bank structure and/or to vegetation growth and its ability to trap and stabilize sediment.                                                                                                                                                                                                                                    	| ~10<sup>-2</sup>     	|
| `b0`         	| **Initial width**. Starting width of a channel                                                                                                                                                                                                                                                                    	| 1&ndash;1000 m                   	|

#### Key input data sets and parameters (FlowDepthDoubleManning)

*This step is used to compute flow depths from a discharge time series, and may be skipped if you already posess a time series of flow depth*

* Discharge time series
* Manning's n (channel)
* Roughness / topogrpahy coefficient (floodplains)
* Depth / topography exponent (floodplains)

### Outputs

This program outputs a time series of channel width, `b(t)`. It organizes this within a Pandas DataFrame that can also be exported using the `write_csv()` function within the `RiverWidth` class.
 
Plots can also be made of just river width (`plotb()`) or of discharge and river width (`plotQb`).
