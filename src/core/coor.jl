## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: Type definition of Coor for coordiates
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================


## ---------------------------------
## define a Super Type for Type 'Coor' and 'Domain'

abstract Location


## ---------------------------------
## define type 'Coor' to hold coordinates

immutable Coor <: Location
    x::Float64
    y::Float64
    time::Float64
end


## !!! outer constructor does not work on current development version of JULIA!!!
## see commented code below !!!

## outer constructor
function Coor(x::Real, y::Real, time::Real)
    Coor(convert(Float64, x),
         convert(Float64, y),
         convert(Float64, time))
end

## set time to 0.0 as default
Coor(x::Real, y::Real) = Coor(x, y, 0.0)

## convert Datetime in milliseconds since REF_TIME
function Coor(x::Real, y::Real, time::DateTime)
    time_msec = time - REF_TIME
    Coor(x, y, time_msec)
end

## !!! Code for development version !!!
## ## set time to 0.0 as default
## Coor(x, y) = Coor(x, y, 0.0)

## ## convert Datetime in milliseconds since REF_TIME
## function Coor(x, y, time::DateTime)
##     time_msec = time - REF_TIME
##     Coor(x, y, time_msec)
## end



## addition and subtraction of coordinates
function +(c1::Coor, c2::Coor)
    Coor(c1.x+c2.x, c1.y+c2.y, c1.time+c2.time)
end

function -(c1::Coor, c2::Coor)
    Coor(c1.x-c2.x, c1.y-c2.y, c1.time-c2.time)
end


## rotate coordinates around 'origin' (in spatial dimensions only)
## angle: rotation angle in rad
function rotate(coor::Coor, origin::Coor, angle::Float64)

    ## short-cut
    angle != 0.0 ? nothing : return(coor)
       
    x_trans = coor.x - origin.x
    y_trans = coor.y - origin.y

    Coor(cos(angle)*x_trans - sin(angle)*y_trans + origin.x,
         sin(angle)*x_trans + cos(angle)*y_trans + origin.y,
         coor.time)
end



## ---------------------------------
## define type 'Domain'

immutable Domain <: Location
    position::Coor
    extend::Coor                    # extend relative to position, before rotation
    angle::Float64                  # angle, rotated around 'position'
end
