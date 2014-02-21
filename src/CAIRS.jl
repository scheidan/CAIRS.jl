## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: Pack everything in a module
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================


## ---------------------------------
## Core module with low-level interface

module CAIRS

## using Distibutions
using Cubature
using Datetime


## -----------
## Core functions

include(joinpath("core", "coor.jl"))                      
include(joinpath("core", "sensor.jl"))
include(joinpath("core", "signal.jl"))
include(joinpath("core", "priorGP.jl"))
include(joinpath("core", "constantes.jl"))
include(joinpath("core", "calib.jl"))
include(joinpath("core", "predict.jl"))

## High-level interface

include(joinpath("interface", "helpers.jl")) #  Helper functions (maybe removed later)
include(joinpath("interface", "signal_handling.jl"))
include(joinpath("interface", "predictions.jl"))            


## -----------
## export

## Types
export Location, Coor, Domain, Sensor, Signal

## core functions
export Gibbs, sample_preditions, sample_realization
export find_near_signals

## interface functions
export add_signal!, remove_signal!, predict
export chains2csv, summary2csv, sensor2csv

## export so that user can supply DateTime objects
export DateTime, datetime


end