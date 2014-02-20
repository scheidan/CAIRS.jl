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
##  Column 1: holds date and time in exactly the following form: "22.11.2013 13:15:30"
##  Column 2: holds the signal values
##
## signal:    Vector with type 'Signal'
## file:      filename as string
## sensor:    a Sensor object
## position:  Coordinates of the sensor, time is ignored
## angle:     angle of sensor, optional
## delim:     delimitation character

function add_signal!(signals::Vector{Signal}, file::String,
                     sensor::Sensor,
                     position::Coor, angle::Float64=0.0;
                     delim=',')

    ## -----------
    ## read file

    data = readdlm(file, delim)

    ## -----------
    ## add signals

    for i in 1:size(data,1)
        ## convert time. Somewhat a hack, no proper parsing yet
        t_string = data[i,1]
        time = datetime(int(t_string[7:10]), int(t_string[4:5]),
                        int(t_string[1:2]), int(t_string[12:13]),
                        int(t_string[15:16]), int(t_string[18:19])) - REF_TIME

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
