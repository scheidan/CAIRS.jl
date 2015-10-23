## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: define type Signal
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

## ---------------------------------
## define type Signal

immutable Signal{T}
    ## signal obtained by sensor
    signal::T

    ## sensor
    sensor::Sensor

    ## coordinates of the sensor
    position::Coor

    ## direction in rad, rotations centre is 'coor'
    angle::Float64

end

## alpha=0 as default
function Signal(signal, sensor, coor)
    Signal(signal, sensor, coor, 0.0)
end

## ---------------------------------
## show

import Base.show
function show(signal::Signal)

    println("Signal:")
    println("- position: $(signal.position)")
    signal.angle != 0.0 ?  println("- angle: ", signal.angle) : nothing
    println("- measured value: $(signal.signal)")
    println("- signal type: $(typeof(signal.signal))")
    println("- sensor:")
    show(signal.sensor, "   ")

end



## ---------------------------------
## Find all signals that are near enough (in time) to
## have an impact on the predicted Locations
## !!! Checks only distance in time !!!
##
## loc_pred:  vector with locations for predictions
## signals:   vector of avaiable Signals
## delta:     max. distance to be considered
##
## returns a vector of signals that are "close" in time to locations

function find_near_signals(loc_pred::Vector, signals::Vector, delta::Real)

    ## Find the extrem time-coordiantes of 'loc_pred'
    tmin = Inf
    tmax = -Inf
    for loc in loc_pred
        if typeof(loc) == Domain
            t1 = loc.position.time
            t2 = (loc.position + loc.extend).time
            tmin = min(tmin, t1, t2)
            tmax = max(tmax, t1, t2)
        else
            tmin = min(tmin, loc.time)
            tmax = max(tmax, loc.time)
        end
    end
    tmin -= delta
    tmax += delta

    ## Find signals that are within (tmin, tmax)
    signals_near = Signal[]
    for sig in signals
        t = [sig.position.time;
             (sig.position + sig.sensor.domain_extent).time;
             [(sig.position + i).time for i in sig.sensor.delta_coor]]

        if any(t .< tmax) && any(t .> tmin)
            push!(signals_near, sig)
        end
    end
    signals_near
end
