## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: tests
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

using Base.Test
using CAIRS

## Coordinates and Domains

@test Coor(1,1,1)==Coor(1.0, 1.0, 1.0)
@test Coor(0,0,0)-Coor(1,1,1)+Coor(1,1,1)==Coor(0,0,0)
@test CAIRS.rotate(CAIRS.rotate(Coor(1,2,3), Coor(0,0,0), 1.0), Coor(0,0,0), -1.0) == Coor(1, 2, 3)

c1 = Coor(1,1,1)
d0 = Domain(c1, Coor(0, 0, 0), 0.0)
d1 = Domain(c1, Coor(-2, 0, 0), 0.0)
d2 = Domain(c1, Coor(2, -4, 0), 0.0)
d3 = Domain(c1, Coor(-2, 4, -8), 0.0)

@test CAIRS.volume(d0) == 0.0
@test CAIRS.volume(d1) == 2.0
@test CAIRS.volume(d2) == 2.0*4.0
@test CAIRS.volume(d3) == 2.0*4.0*8.0

## sensors

foo(x) = x
@test Sensor(foo).delta_coor == [Coor(0,0,0)] # point measurement
@test Sensor(foo).domain_extent == Coor(0,0,0)

@test Sensor(foo, Coor(1,1,1)).domain_extent == Coor(1,1,1) # not offset, only integration domain
@test Sensor(foo, Coor(1,1,1)).delta_coor == Coor[]


## Run example script

module examtest
include(joinpath(Pkg.dir("CAIRS"), "example", "Example_highlevel_interface.jl"))
end
