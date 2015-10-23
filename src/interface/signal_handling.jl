## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: high-level interface to import signals
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================



## ---------------------------------
## add signals of a file to a signal vector
##
## The file must contain two columns:
##  Column 1: holds date and time
##  Column 2: holds the signal values
##
## signal:      vector with type 'Signal'
## file:        filename as string
## sensor:      a Sensor object
## position:    coordinates of the sensor, time is ignored
## angle:       angle of sensor, optional
## delim:       delimitation character of file
## date_format: format string to parse dates and time

function add_signal!(signals::Vector{Signal}, file::AbstractString,
                     sensor::Sensor,
                     position::Coor, angle::Float64=0.0;
                     delim::Char=',',
                     date_format::AbstractString="d.m.yyyy HH:MM:SS")

    ## -----------
    ## read file

    data = readdlm(file, delim)

    ## -----------
    ## add signals

    for i in 1:size(data,1)
        time = DateTime(data[i,1], date_format)
        new_signal = Signal(data[i, 2], sensor, Coor(position.x, position.y, time), angle)
        push!(signals, new_signal)

    end

end



## -----------
## remove signals from the signal vector

## signal:    Vector with type 'Signal'
## sensor:    a Sensor object
## position:  Coordinates of the sensor, time is ignored
## angle:     angle of sensor, optional


## remove signals measured with "sensor"
function remove_signal!(signals::Vector{Signal},
                        sensor::Sensor)

    filter!(x -> x.sensor!=sensor, signals)
end

## remove signals with this position and angle
function remove_signal!(signals::Vector{Signal},
                        position::Coor,
                        angle::Float64=0.0)

    f_filter(x) = !(x.position.x==position.x && x.position.y==position.y && x.angle==angle)
    filter!(f_filter, signals)
end
