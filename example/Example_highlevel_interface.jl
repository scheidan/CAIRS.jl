## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: minimal example of the high-level interface
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================


using CAIRS
using Datetime

## ---------------------------------
## define sensors

## -- non-linear continous sensor
function log_p_gauge(S::Float64, R::Vector)

    mu = 0.1+R[1]^2.0
    sigma = 0.005

    ## log of normal density
    -(S-mu)^2.0/(2.0*sigma)
end

sensor_gauge = Sensor(log_p_gauge, [Coor(0.0, 0.0, -5*60*1000.0)]) # integrates 5 minutes in time

## -- microwave link
function log_p_MWL(S::Float64, I::Float64)

    R_mean = I/6.0
    sigma = 0.005

    ## log of normal density
    -(S-R_mean)^2.0/(2.0*sigma)
end

sensor_MWL = Sensor(log_p_MWL, Coor(6.0, 0.0, 0.0)) # integrates along a path


## ---------------------------------
## Import signals from files

## The file must contain two columns:
##  Column 1: holds date and time in exactly the following form: "22.11.2013 13:15:30"
##  Column 2: holds the signal values


sig = Signal[]                          # create an empty array

add_signal!(sig, joinpath(Pkg.dir("CAIRS"), "example/data/Sensor1.csv"), sensor_gauge, Coor(5, 6), delim=',')
add_signal!(sig, joinpath(Pkg.dir("CAIRS"), "example/data/Sensor2.csv"), sensor_MWL, Coor(4.2, 2), 0.9)

show(sig[end])

## -- remove signals if necessary
## remove_signal!(sig, sensor_gauge)
## remove_signal!(sig, Coor(10, 20, -12))


## write sensor position in a file (useful for plotting)
sensor2csv(sig, "sensor_coor.csv")


## -----------
## Define location for predictions

## create a simple grid (iregular predictions are possible too)
nn = 21
loc_pred = [Coor(i, j, time)
            for i=linspace(0, 10, nn), j=linspace(0, 10, nn),
            time=datetime(2013, 11, 22, 13, 16, 00) : minute(1): datetime(2013, 11, 22, 13, 20, 00) ]


## -----------
## find near signals

R_pred = predict(loc_pred, sig,
                 n_sample_calib = 30000,
                 burn_in = 1200,
                 n_sample_pred = 6000,
                 delta = 60*1000)

## compute  and plot (with R)
summary2csv(R_pred, "rain_field.csv")

## call R to make a plot
## run(`Rscript $(joinpath(Pkg.dir("CAIRS"), "R", "compute_rain_map.r")) rain_field.csv sensor_coor.csv out.pdf`)
)
