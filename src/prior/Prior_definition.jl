## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: mean and covariance functions to define prior
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

## Functions that return a mean or covariance function for arguments
## of type "Coor".


## ---------------------------------
## mean functions

## -----------
## constant

function mean_constant(;mean::Real=0.0)

    println("- Prior has a constant mean of $mean.")

    function f_mu(c::Coor)
        mean
    end

    return f_mu
end




## ---------------------------------
## covariance (kernel) function

## -----------
## separable gamma-exponential

## Equation (4.18) in Rasmussen, C. E. and Williams, C. K. I. (2006) Gaussian processes
## for machine learning, Massachusett, MIT Press.

function cov_exponential(;sigma::Real=10.0,       # standard deviation of GP
                         l_spatial::Real=3000.0, # spatial correlation length
                         l_temporal::Real=600.0*1000, # temporal correlation length [milliseconds]
                         gamma::Real=1.0)        # exponent for smoothness in [0, 2]


    (gamma<0) | (gamma>2) ? error("'gamma' must be in [0, 2]!") : nothing

    println("- Prior has a separable gamma-exponential covariance function with:")
    println("   standard deviation: $sigma")
    println("   spatial correlation length: $l_spatial")
    println("   temporal correlation length [ms]: $l_temporal")
    println("   gamma: $gamma")


    var = sigma^2

    function f_cov(c1::Coor, c2::Coor)

        dist_spatial = sqrt((c1.x-c2.x)^2 + (c1.y-c2.y)^2)
        dist_temoral = abs(c1.time-c2.time)

        var * exp( -(dist_spatial/l_spatial)^gamma - (dist_temoral/l_temporal)^gamma )
    end

    return f_cov
end
