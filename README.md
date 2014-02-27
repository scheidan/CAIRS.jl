# CAIRS - Continuous Assimilation of Integrating Rain Sensors


![rain map](https://raw.github.com/scheidan/CAIRS/master/Images%20for%20Readme/Header.png)

[![Build Status](https://travis-ci.org/scheidan/CAIRS.jl.png)](https://travis-ci.org/scheidan/CAIRS.jl)

_CAIRS_ is a framework to reconstruct rain fields by assimilating
signals of fundamentally different rain sensors .

In particular, the *integration characteristics* of sensors are
explicitly considered.  For example, non-recording standard rain gauges
integrate over time and deliver information such as the daily rainfall
sums. The rain-induced attenuation of micro wave links (MWL) can be
used to measure the path-integrated intensities - an example of
spatial integration.

Sensor signals with different scales (e.g. continuous, binary) can be
assimilated. Furthermore, _CAIRS_ is formulated continuously in time
and space. This is helpful because it enables a natural consideration
of signals with irregular time-intervals.

For more information see [Scheidegger and Rieckermann (2014)](#publication).

Note, _CAIRS_ is still in development and the interface may change.



# Installation

_CAIRS_ is a [Julia](http://julialang.org/) package. The first step is to download and install
Julia (http://julialang.org/downloads/).

_CAIRS_ can then be installed easely with the Julia command `Pkg.clone()`:

```Julia
Pkg.clone("git://github.com/scheidan/CAIRS.jl.git")
```

After that, _CAIRS_ behaves like a normal package. For example, it can
be updated with `Pkg.update()`.


# Example

First, the package _CAIRS_ must be loaded. For convinience, it is also
recommended to load the packages `Datetime` and `Distributions`:

```Julia
using CAIRS
using Datetime
using Distributions
```


### Sensor definition

Every sensor must be defined. In the simplest case a sensor measures
the rain intensity at a point. Then simply (the logarithm of) the signal
distribution must be defined:

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

### Signal import

The next step is to import the signals. Every signal must have an
attached attached sensor. Signals can be constructed with the function
`Signal` or more conveniently with `add_signal()`.

Currently `add_signal()' expected that the signals of every sensor are
stored in a separate file. The file must contain two columns:
- Column 1: date and time in *exactly* the following form: "22.11.2013 13:15:30"
- Column 2: signal values

```julia
## path to example data that come with the CAIRS package
path1 = joinpath(Pkg.dir("CAIRS"), "example", "data", "Sensor1.csv")
path2 = joinpath(Pkg.dir("CAIRS"), "example", "data", "Sensor2.csv")

sig = Signal[]                          # create an empty array for Signals

add_signal!(sig,                        # add signal to vector 'sig'
            path1,                      # file name
            sensor_gauge,               # gauge sensor, defined above
            Coor(5, 6)                  # coordinate of the sensor
            )                           # optional argument: delim=','

add_signal!(sig, path2,
            sensor_MWL,                 # MWL sensor, defined above
            Coor(4.2, 2),               # coordinate of one end point of the sensor
            0.9)                        # rotation around the point defined above in [rad]
```

Information about a signal can be printed with `show`, e.g. `show(sig[1])`.

The sensor positions can be written to a file with:
```Julia
sensor2csv(sig, "sensor_coor.csv")
```
This is useful for plotting.

### Definition of prediction points

The location for which a prediction of the rain intesity is desired must be
defined as an `Array` or `Vector` of coordinates. Coordinates are
defined with `Coor(x, y, time)`. Time can be a number or a `datetime`
object.
```Julia
### create a simple grid
nn = 20
loc_pred = [Coor(i, j, time)
            for i=linspace(0, 10, nn), j=linspace(0, 10, nn),
            time=datetime(2013, 11, 22, 13, 15, 00) : minute(1): datetime(2013, 11, 22, 13, 20, 00) ]
```
This produced a regular grid, but the point could also be irregularly distributed.

### Assimilation
The assimilation of the signals and the computation of the predictions are done with `predict`.
```Julia
R_pred = predict(loc_pred,               # vector or array with locations for predictions
                 sig,                    # vector of signals
                 n_sample_calib = 20000, # number of iterations of the Gibbs sampler
                 burn_in = 5000,         # number of removed samples (and length of adaptation)
                 n_sample_pred = 6000,   # number of samples for predictions
                 delta = 90*1000)        # consider all signals that are not further away than
                                         # time 'delta' from prediction points [ms]
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
Note, here is is assumed that `Rscript` is in PATH.



# Publication
<a name="Publication"></a>

Scheidegger, A. and Rieckermann, J. (2014) "Bayesian assimilation of
rainfall sensors with fundamentally different integration
characteristics" in WRaH Proceedings, Washington, DC.
