## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: Type definition of Sensors
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

## ---------------------------------
## define Sensor

struct Sensor
    ## function that computes log p(S| [R_1 ,... R_n], I)
    ## must take a signal S and an array of *not transformed* rain
    log_p::Function

    ## 'offset' in the coordinates of all [R_1 ,... R_n]
    delta_coor::Vector{Coor}

    ## cube of integration, relative to sensor position
    domain_extent::Coor
end

## constructor assumes no integration if not specified
function Sensor(log_p::Function, delta_coor::Vector{Coor})
    coor = Coor(0.0, 0.0, 0.0)
    Sensor(log_p, delta_coor, coor)
end

## constructor assumes no offset if not specified
function Sensor(log_p::Function, int_domain::Coor)
    Sensor(log_p, Coor[], int_domain)
end

## constructor assumes no offset and no integration if not specified
function Sensor(log_p::Function)
    coor = Coor(0.0, 0.0, 0.0)
    Sensor(log_p, [coor], coor)
end


## ---------------------------------
## show

import Base.show
function show(sensor::Sensor, offset=""::String)

    length(offset)==0 ? println("Sensor:") : nothing
    println(offset,"- log likelihood: $(sensor.log_p)")

    if size(sensor.delta_coor, 1) > 0
        print(offset, "- signal is condition on:\n")
        for coor in sensor.delta_coor
            println(offset, "   ", coor)
        end
    end

    if sensor.domain_extent != Coor(0.0, 0.0, 0.0)
        println(offset, "- integration domain: $(sensor.domain_extent)")
    else
        println(offset, "- no integration")
    end

end
