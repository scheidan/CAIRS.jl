# CAIRS - Continuous Assimilation of Integrating Rain Sensors


![rain map](https://raw.github.com/scheidan/CAIRS/master/Images%20for%20Readme/Header.png)

Linux, OS X:
[![Build Status](https://travis-ci.org/scheidan/CAIRS.jl.svg?branch=master)](https://travis-ci.org/scheidan/CAIRS.jl)
  Windows:
[![Build status](https://ci.appveyor.com/api/projects/status/9cuvjuek83wut0ju/branch/master?svg=true)](https://ci.appveyor.com/project/scheidan/cairs-jl/branch/master)
[![Coverage Status](https://img.shields.io/coveralls/scheidan/CAIRS.jl.svg)](https://coveralls.io/r/scheidan/CAIRS.jl?branch=master)

_CAIRS_ is a framework to reconstruct rain fields by assimilating
signals of fundamentally different rain sensors.

In particular, the *integration characteristics* of sensors are
explicitly considered.  For example, non-recording standard rain gauges
integrate over time and deliver information such as the daily rainfall
sums. The rain-induced attenuation of micro wave links (MWL) can be
used to measure the path-integrated intensities - an example of a sensor with
spatial integration.

Sensor signals with different scales (e.g. continuous, binary) can be
assimilated. Furthermore, _CAIRS_ is formulated continuously in time
and space. This is helpful because it enables a natural consideration
of signals with irregular time-intervals.

The mathematical model is described in [Scheidegger and Rieckermann (2014)](#publication).
The basic functionality and application is explained in [this tutorial](https://github.com/scheidan/CAIRS-Tutorial).

Note, _CAIRS_ is still under development and the interface may change.



# Installation

_CAIRS_ is a [Julia](http://julialang.org/) package. The first step is to download and install
[Julia](http://julialang.org/downloads/) version 0.5 or newer.

_CAIRS_ can then be installed with the Julia command `Pkg.clone()`:

```Julia
Pkg.clone("git://github.com/scheidan/CAIRS.jl.git")
```

After that, _CAIRS_ behaves like a normal package. For example, it can
be updated with `Pkg.update()`.


# Example

First, the package _CAIRS_ must be loaded. For convinience, it is also
recommended to load the packages `Dates` and `Distributions`:

```Julia
using CAIRS
using Distributions
using Base.Dates
```


### Sensor definition

Every sensor must be characterized. In the simplest case a sensor measures
the rain intensity at a point. In this case (the logarithm of) the signal
distribution must be defined conditioned on the intesity at this coordinate:

```Julia
function log_p_gauge(S::Float64, R::Vector) # non-linear continuous rain gauge

	mu = 0.1+R[1]^2.0    # Note, the signal and can be non-linearly
						 # related to the rain intensity.
	sigma = 0.005

	 ## log of normal density, p(S|R)
	logpdf(Normal(mu, sigma), S)   # doesn't have to be normal
end

sensor_gauge = Sensor(log_p_gauge)
```

For integrating sensors, also the integration domain must be
specified. For example, a micro wave link (MWL) with length `6` may be
defined as:

```julia
function log_p_MWL(S::Float64, I::Float64)

	R_mean = I/6.0
	sigma = 0.1

	## log of normal density, p(S|I)
	logpdf(Normal(R_mean, sigma), S)
end

sensor_MWL = Sensor(log_p_MWL, Coor(6, 0, 0)) # integrates along a path of length 6
```

### Prior definition

The prior of the rain field is modeled as Gaussian process (GP). A GP
is described by a mean and a covariance function.

This functions can be specified by the user. The mean function returns
the prior mean of the rain intensity at a given coordinate. It must
take a single argument of type `Coor`. The covariance function must
return the covariance of the rain intensities at two given point, given
by two arguments of type `Coor`. Note, it is not checked if the
provided function is a valid covariance function!

However, helpers to construct valid functions are provided. The functions
`mean_constant()` and `cov_exponential()` create a simple constant
mean, and a separable gamma-exponential covariance function. Only the
parameters must be provided:

```Julia
mean_GP = mean_constant(mean=2.0)

cov_GP = cov_exponential(sigma=10.0,           # standard deviation of GP
						 l_spatial=1.5,        # spatial correlation length
						 l_temporal=Minute(1), # temporal correlation length
						 gamma=1.0)            # exponent for smoothness in [0, 2]
```
Other types of covariance functions will be added in future.

### Signal import

The next step is to import the signals. Every signal must have an
attached sensor. Signals can be constructed with the function
`Signal` or more conveniently with `add_signal()`.

Currently `add_signal()` expected that the signals of every sensor are
stored in a separate file. The file must contain two columns:

 - Column 1: date and time
 - Column 2: signal values

```julia
## path to example data that come with the CAIRS package
path1 = joinpath(Pkg.dir("CAIRS"), "example", "data", "Sensor1.csv")
path2 = joinpath(Pkg.dir("CAIRS"), "example", "data", "Sensor2.csv")

sig = Signal[]                          # create an empty array for Signals

add_signal!(sig,                        # add signal to vector 'sig'
			path1,                      # file name
			sensor_gauge,               # sensor
			Coor(5, 6),                 # coordinate of the sensor
			date_format="d.m.yyyy HH:MM:SS",
			delim=',')                  # delimitation character


add_signal!(sig, path2,
			sensor_MWL,                 # MWL link
			Coor(4.2, 2),               # coordinate of one end point of the sensor
			0.9,                        # rotation around the point defined above in [rad]
			date_format="d.m.yyyy HH:MM:SS",
			delim=',')
```

Information about a signal can be printed with `show`, e.g. `show(sig[1])`.

Writing the sensor positions in a file is useful for plotting:
```Julia
sensor2csv(sig, "sensor_coor.csv")
```

### Definition of prediction points

The location for which a prediction of the rain intesity is desired must be
defined as an `Array` or `Vector` of coordinates. Coordinates are
defined with `Coor(x, y, time)`. Time can be a number or a `DateTime`
object.
```Julia
### create a simple grid
nn = 20
loc_pred = [Coor(i, j, time)
			for i=linspace(0, 10, nn), j=linspace(0, 10, nn),
			time=DateTime(2013, 11, 22, 13, 15, 00) : Minute(1): DateTime(2013, 11, 22, 13, 20, 00) ]
```
This produced a regular grid, but the point could also be irregularly distributed. Also, not only predictions for coordinates but also for intesities integrated over a domain can be made. Domains are defined by the function `Domain`.

### Assimilation
The assimilation of the signals and the computation of the predictions are done with `predict`.
```Julia
R_pred = predict(loc_pred,               # vector or array with locations for predictions
				 sig,                    # vector of signals
				 mean_GP,                # mean function of prior
				 cov_GP,                 # covariance function of prior
				 n_sample_calib = 20000, # number of iterations of the Gibbs sampler
				 burn_in = 5000,         # number of removed samples (and length of adaptation)
				 n_sample_pred = 6000,   # number of samples for predictions
				 delta = Second(90))     # consider all signals within time 'delta'
										 # from prediction points
```

Write a summary of the samples in a file that is used for visualization:
```julia
summary2csv(R_pred, "rain_field.csv")
```

### Visualization with R
One possibility to visualize the result is to use [R](http://www.r-project.org/). A simple
R-script to produce rain maps comes with _CAIRE_. It requires that R and
the R-libraries `lattice`, `latticeExtra` and `tripack` are installed.
```Julia
pathRscript = joinpath(Pkg.dir("CAIRS"), "R", "compute_rain_map.r")
run(`Rscript $pathRscript  rain_field.csv sensor_coor.csv out.pdf`)
```
Note, here it is assumed that `Rscript` is in PATH.



# Publication
<a name="Publication"></a>

Scheidegger, A. and Rieckermann, J. (2014) "Bayesian assimilation of
rainfall sensors with fundamentally different integration
characteristics" in WRaH Proceedings, Washington, DC.
