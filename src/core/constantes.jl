## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: Contains all constant variables
##
## File: import.jl
## Path: c:/Users/scheidan/Dropbox/Eawag/Rainfall assimilation/Julia/CRA/
##
## November 14, 2013 -- Andreas Scheidegger
##
## andreas.scheidegger@eawag.ch
## =======================================================

## -----------
## parameters of covariance function

const var = 10*10                       # variance of GP
const l_spatial = 1500.0                # spatial correlation length [m]
const l_temporal = 50*60*1000.0        # temporal correlation length [milli seconds]
const p = 1.0                           # exponent for smoothness in [0, 2]
const alpha = 0.5                       # weight between spatial and temporal cor.

## -----------
## misc

const REF_TIME = datetime(1984, 10, 20) # reference point in time. Count milliseconds from this date.
