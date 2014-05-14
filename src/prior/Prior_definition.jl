## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: mean and covariance functions to define prior
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================



## ------------------------------------------------------------------
## transformation
##
## if X~N(0,1), g(X) should have a realistic rain distibution

## exponential with offset
function trans2real(R_trans::Float64)
    k = 0.8
    offset = 1.5
    R_trans <= offset ? 0 : exp(k*(R_trans-offset)) - 1.0
end

## ## power law with limitation and offset
## function trans2real(R_trans::Float64)
##     k = 2.0
##     offset = 1.8
##     a = 3.5
##     if R_trans <= offset
##         return(0.0)
##     else
##         if(R_trans < a)
##             return( (R_trans-offset)^k )
##         else
##             slope = k*(a-offset)^(k-1)
##             return( (R_trans-a)*slope + (a-offset)^k )
##         end
##     end
## end

## ## no transformation
## function trans2real(R_trans::Float64)
##    R_trans
## end

## vectorized version
function trans2real(R_trans::Vector{Float64})
    n = size(R_trans,1)
    R_real = zeros(n)
    for i in 1:n
        R_real[i] = trans2real(R_trans[i])
    end
    return(R_real)
end

## ------------------------------------------------------------------
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
## covariance (kernel) functions

## -----------
## separable gamma-exponential

## Equation (4.18) in Rasmussen, C. E. and Williams, C. K. I. (2006) Gaussian processes
## for machine learning, Massachusett, MIT Press.

function cov_exponential(;sigma::Real=10.0,       # standard deviation of GP
                         l_spatial::Real=3000.0, # spatial correlation length
                         l_temporal=600.0*1000, # temporal correlation length [milliseconds]
                         gamma::Real=1.0)        # exponent for smoothness in [0, 2]

    (gamma<0) | (gamma>2) ? error("'gamma' must be in [0, 2]!") : nothing

    ## convert TimePeriod objects in milliseconds
    if typeof(l_temporal) <: Datetime.TimePeriod
        temp = l_temporal - second(0)
        l_temporal = 1000*convert(Real, temp)
    end

    println("- Prior has a separable gamma-exponential covariance function with:")
    println("   standard deviation: $sigma")
    println("   spatial correlation length: $l_spatial")
    println("   temporal correlation length [s]: $(l_temporal/1000)")
    println("   gamma: $gamma")

    var = sigma^2

    function f_cov(c1::Coor, c2::Coor)
        dist_spatial = sqrt((c1.x-c2.x)^2 + (c1.y-c2.y)^2)
        dist_temoral = abs(c1.time-c2.time)
        var * exp( -(dist_spatial/l_spatial)^gamma - (dist_temoral/l_temporal)^gamma )
    end

    return f_cov
end


## -----------
## local spherical (compact support)

## Modern Applied Statistics with S (2002) New York, Springer, page 428.
## Gneiting, T. (2002) Compactly Supported Correlation Functions. Journal
##    of Multivariate Analysis, 83(2), 493-508.

function cov_sphere(;sigma::Real=10.0,       # standard deviation of GP
                         l_spatial::Real=3000.0, # spatial correlation length
                         l_temporal=600.0*1000) # temporal correlation length [milliseconds]

    ## convert TimePeriod objects in milliseconds
    if typeof(l_temporal) <: Datetime.TimePeriod
        temp = l_temporal - second(0)
        l_temporal = 1000*convert(Real, temp)
    end

    println("- Prior has a spherical covariance function with:")
    println("   standard deviation: $sigma")
    println("   spatial radius: $l_spatial")
    println("   temporal radius [s]: $(l_temporal/1000)")

    var = sigma^2
    function f_cov(c1::Coor, c2::Coor)
        r = sqrt( ((c1.x-c2.x)/l_spatial)^2 + ((c1.y-c2.y)/l_spatial)^2 + ((c1.time-c2.time)/l_temporal)^2 )
        if r>1.0
            cov = 0.0
        else
            cov = var*(1-3/2*r + 1/2*r^3)
        end
        cov
    end

    return f_cov
end
