## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: tests
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

using Base.Test
using CAIRS


## Coordinates

@test Coor(1,1,1)==Coor(1.0, 1.0, 1.0)
@test Coor(0,0,0)-Coor(1,1,1)+Coor(1,1,1)==Coor(0,0,0)
@test CAIRS.rotate(CAIRS.rotate(Coor(1,2,3), Coor(0,0,0), 1.0), Coor(0,0,0), -1.0) == Coor(1, 2, 3)


## sensors

foo(x) = x
@test Sensor(foo).delta_coor == [Coor(0,0,0)] # point measurement
@test Sensor(foo).domain_extent == Coor(0,0,0)

@test Sensor(foo, Coor(1,1,1)).domain_extent == Coor(1,1,1) # not offset, only integration domain
@test Sensor(foo, Coor(1,1,1)).delta_coor == Coor[]


## Signals
