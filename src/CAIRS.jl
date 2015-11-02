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
using Base.Dates

import PDMats

import Base.+
import Base.-

## -----------
## Define constantes

const REF_TIME = Dates.DateTime(1984, 10, 20) # reference point in time. Count milliseconds from this date.

## -----------
## Core functions

include(joinpath("core", "coor.jl"))
include(joinpath("core", "sensor.jl"))
include(joinpath("core", "signal.jl"))
include(joinpath("core", "GaussianProcess.jl"))
include(joinpath("core", "calib.jl"))
include(joinpath("core", "predict.jl"))

## High-level interface

include(joinpath("interface", "helpers.jl")) #  Helper functions (maybe removed later)
include(joinpath("interface", "signal_handling.jl"))
include(joinpath("interface", "predictions.jl"))


## Prior
include(joinpath("prior", "Prior_definition.jl"))

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

## prior
export mean_constant, cov_exponential, cov_sphere


end
