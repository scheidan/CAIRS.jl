## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: Type definition of Sensors
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

## ---------------------------------
## define type Sensor

immutable Sensor
    ## function that computes log p(S| [R_1 ,... R_n], I)
    ## must take a signal S and an array of *not transformed* rain
    log_p::Function

    ## 'offset' in the coordinates of all [R_1 ,... R_n]
    delta_coor::Vector{Coor}

    ## cube of integration, relative to sensor position
    domain_extent::Coor

    ## function applied on the rain before integration: \int f_int(R(c)) dc
    f_int::Function
end

## constructor assumes no integration if not specified
function Sensor(log_p::Function, delta_coor::Vector{Coor})
    coor = Coor(0.0, 0.0, 0.0)
    identity(x) = x
    Sensor(log_p, delta_coor, coor, identity)
end

## constructor assumes no integration if not specified
function Sensor(log_p::Function, delta_coor::Vector{Coor})
    coor = Coor(0.0, 0.0, 0.0)
    identity(x) = x
    Sensor(log_p, delta_coor, coor, identity)
end

## constructor assumes no offset if not specified
function Sensor(log_p::Function, int_domain::Coor)
        identity(x) = x
    Sensor(log_p, Coor[], int_domain, identity)
end

## constructor assumes no offset if not specified
function Sensor(log_p::Function, int_domain::Coor, f_int::Function)
    Sensor(log_p, Coor[], int_domain, f_int)
end

## constructor assumes no offset if not specified
function Sensor(log_p::Function, delta_coor::Vector{Coor}, int_domain::Coor)
    identity(x) = x
    Sensor(log_p, delta_coor, int_domain, identity)
end

## constructor assumes no offset and no integration if not specified
function Sensor(log_p::Function)
    coor = Coor(0.0, 0.0, 0.0)
    identity(x) = x
    Sensor(log_p, [coor], coor, identity)
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
        println(offset, "- integration function: $(sensor.f_int)")
    else
        println(offset, "- no integration")
    end

end
