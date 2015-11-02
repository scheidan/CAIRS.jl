## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: minimal example of the high-level interface
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================


using CAIRS
using Base.Dates
using Distributions

## ---------------------------------
## Define sensors

## -- non-linear continuous rain gauge
function log_p_gauge(S::Float64, R::Vector)

    mu = 0.1+R[1]^2.0
    sigma = 0.005

    ## log of normal density, p(S|R)
    logpdf(Normal(mu, sigma), S)
end

sensor_gauge = Sensor(log_p_gauge) # no integration

## -- microwave link
function log_p_MWL(S::Float64, I::Float64)

    R_mean = I/6.0
    sigma = 0.005

    ## log of normal density, p(S|I)
    logpdf(Normal(R_mean, sigma), S)
end

sensor_MWL = Sensor(log_p_MWL, Coor(6, 0, 0)) # integrates along a path of length 6


## ---------------------------------
## Define prior

mean_GP = mean_constant(mean=2.0)

cov_GP = cov_exponential(sigma=10.0,           # standard deviation of GP
                         l_spatial=1.5,        # spatial correlation length
                         l_temporal=Minute(1), # temporal correlation length
                         gamma=1.0)           # exponent for smoothness in [0, 2]


## ---------------------------------
## Import signals from files

## Signals of every sensor must be in a separate file.
## The file must contain two columns:
##  Column 1: holds date and time
##  Column 2: holds the signal values


sig = Signal[]                          # create an empty array

## path to example data
path1 = joinpath(Pkg.dir("CAIRS"), "example", "data", "Sensor1.csv")
path2 = joinpath(Pkg.dir("CAIRS"), "example", "data", "Sensor2.csv")

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


## show details of a signal object
show(sig[end])

## -- remove signals if necessary
## remove_signal!(sig, sensor_gauge)
## remove_signal!(sig, Coor(10, 20, -12))


## write sensor position in a file (useful for plotting)
sensor2csv(sig, "sensor_coor.csv")


## -----------
## Define location for predictions

## create a simple grid (irregular predictions are possible too)
nn = 20
loc_pred = [Coor(i, j, time)
            for i=linspace(0, 10, nn), j=linspace(0, 10, nn),
            time=DateTime(2013, 11, 22, 13, 15, 00) : Minute(1): DateTime(2013, 11, 22, 13, 20, 00) ]


## -----------
## Compute predictions

@time R_pred = predict(loc_pred,               # vector or array with locations for predictions
                 sig,                    # vector of signals
                 mean_GP,                # mean function of prior
                 cov_GP,                 # covariance function of prior
                 n_sample_calib = 20000, # number of iterations of the Gibbs sampler
                 burn_in = 5000,         # number of removed samples (and length of adaptation)
                 n_sample_pred = 6000,   # number of samples for predictions
                 delta = Second(90))    # consider all signals within time 'delta'
                                         #  from prediction points


## compute summary of samples and save in file
summary2csv(R_pred, "rain_field.csv")


## -----------
## Visualize the result with R

## One possibility to visualize the result is to use R. A simple
## R-script comes with the package.  R (http://www.r-project.org/) and
## the R-libraries 'lattice', 'latticeExtra' and 'tripack' must be
## installed.

pathRscript = joinpath(Pkg.dir("CAIRS"), "R", "compute_rain_map.r")
run(`Rscript $pathRscript  rain_field.csv sensor_coor.csv out.pdf`)

## here it is assumed that 'Rscript' is in PATH.
