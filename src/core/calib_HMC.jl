## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: calibration of rainfield with Hamiltonian Monte Carlo
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================




## --- structure to implement ? ---
##
## 1) loop through all Signals:
##     i) approximate Domains with coordinates,
##    ii) add new coordinates to in a Vector 'r_calib' if necessary
##   iii) build a dict: Signal => [index of relevant coordinates in 'r_calib']
##
## 2) function to compute log joint prob p(S, R), argument "r_calib":
##      sum( signal.sensor.log_p(S, r_calib[Dict[signal]]) ) + log.prior(r_calib)
##
## 3) construct gradient function of log  p(S, R)
##     - augment each Sensor object with 'log_p_derivative()'
##       Then gradient:
##            r1:  log_p_derivative_s1([r1]) + log_p_derivative_s3([r1, r2]) + log_p_derivative_prior([r1, ...rn])
##       Would need a Dict: r -> [sig1, sig2]
##
##      Maybe with a macro?
##
##      AD-packages:
##       Pkg. ReverseDiffSource.jl, works with most distribution (in combination with MCMC.jl)
##       ForwardDiff.jl
##
##
## 4) Use MCMC package for inference


## ---------------------------------
## Prepares MCMC
## Returns a tuple with:
## - Discretize integration domain -> coors
## - Produce a vector with the coordinates where rain intensity will
##    be sampled -> dic_delta_coor_index
## - Produce two dictionaries with the indices that point to the right
##    rain intensities for each signal -> dic_domain_coor_index

## signals:    vector of Signals
## n_approx:   number of Coor along one dimension to represent a domain

function setupMCMC{T<:Signal}(signals::Vector{T}, n_approx::Integer)

    coors = Coor[]                      # array holds all coordinates that are 'measured'
    dic_delta_coor_index = [S => Int64[] for S in signals]  # Signal => [indices of coors directly measured]
    dic_domain_coor_index = [S => Int64[] for S in signals]  # Signal => [indices of coors for approx. domains]

    ## loop over all Signals
    for S in signals
        ## -- point measurements
        coors_delta = [c + S.position for c in S.sensor.delta_coor]

        ## construct index for direct measurements
        for coor in coors_delta
            if in(coor, coors) # Coor is already in coors
                append!(dic_delta_coor_index[S], findin(coors, [coor])) # add position to index
            else
                append!(coors, [coor]) # add coor
                append!(dic_delta_coor_index[S], [size(coors, 1)] )# add position to index
            end
        end

        ## -- integrated measurements
        ## approximate domains with points
        if S.sensor.domain_extent != Coor(0.0, 0.0, 0.0)

            ## get points for approximation
            coors_approx = domain2points(S, n_approx)

            ## construct index for domains
            for coor in coors_approx
                if in(coor, coors) # coor is already in coors
                    append!(dic_domain_coor_index[S], findin(coors, [coor])) # add position to index
                else
                    append!(coors, [coor]) # add Coor
                    append!(dic_domain_coor_index[S], [size(coors, 1)] )# add position to index
                end
            end
        end

    end

    return(coors, dic_delta_coor_index, dic_domain_coor_index)
end


## ---------------------------------
## defines points to approximate a domain
##
## ??? cleverer approximation scheme ???
## e.g : http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html
##
## S:        Signal with a domain to approximate

## !!! 'cube' in log_likeli() must be adapted if this function is changed !!!

function domain2points(S::Signal, n_approx::Integer)

    ## distribute points along one dimension
    dist_points_1d(extent::Real) = unique(linspace(extent/(2*n_approx), extent*(1-1/(2*n_approx)), n_approx))

    ## rotate and shift in right position
    coors = [rotate(Coor(x, y, time), Coor(0.0, 0.0, 0.0), S.angle) + S.position
             for x in dist_points_1d(S.sensor.domain_extent.x),
             y in dist_points_1d(S.sensor.domain_extent.y),
             time in dist_points_1d(S.sensor.domain_extent.time)]

    coors = reshape(coors, length(coors)) # reshape to vector
    return(coors)
end


## ---------------------------------
## computes log likelihood p(S|R)

## R:        vector with rain intensities corresponding to coors
## signals:  vector of 'Signal' object
## dic_delta_coor_index:     dictionaries as created by setupMCMC()
## dic_domain_coor_index:    dictionaries as created by setupMCMC()
## n_approx:   number of Coor along one dimension to represent a domain

function log_likeli{T<:Signal}(R::Vector{Float64}, signals::Vector{T},
                               dic_delta_coor_index, dic_domain_coor_index,
                               n_approx::Integer)

    ## compute Volume that is represented by one coordinate
    function cube(signal::Signal, n_approx::Integer)
        (signal.sensor.domain_extent.x/n_approx)^(0!=signal.sensor.domain_extent.x) *
        (signal.sensor.domain_extent.y/n_approx)^(0.0!=signal.sensor.domain_extent.y) *
        (signal.sensor.domain_extent.time/n_approx)^(0.0!=signal.sensor.domain_extent.time)
    end

    ll = 0
    for S in signals

        if S.sensor.domain_extent == Coor(0.0, 0.0, 0.0) # no domain
            ll += S.sensor.log_p(S.signal, trans2real(R[dic_delta_coor_index[S]]))
        end
        if size(S.sensor.delta_coor,1) ==  0 # no coordiantes

            int_cube = cube(S, n_approx)
            I = sum(S.sensor.f_int(trans2real(R[dic_domain_coor_index[S]])))*int_cube

            ll += S.sensor.log_p(S.signal, I)
        end
        if (S.sensor.domain_extent != Coor(0.0, 0.0, 0.0)) && (size(S.sensor.delta_coor,1) > 0)

            int_cube = cube(S, n_approx)
            I = sum(S.sensor.f_int(trans2real(R[dic_domain_coor_index[S]])))*int_cube

            ll += S.sensor.log_p(S.signal, trans2real(R[dic_delta_coor_index[S]]), I)
        end

    end
    return(ll)
end



## ---------------------------------
## computes log of p(S, R)
function log_joint{T<:Signal}(R::Vector{Float64}, signals::Vector{T},
                              dic_delta_coor_index, dic_domain_coor_index,
                              mu, Sigma_inv,
                              n_approx)

    log_likeli(R, signals, dic_delta_coor_index, dic_domain_coor_index, n_approx) + log_p_prior(R, mu, Sigma_inv)
end

## ---------------------------------
## run sampler based on package MCMC
##
## signals:      vector of 'Signals'
## prior_mean:   mean function of Prior, f(c::Coor)
## prior_cov:    covariance function of Prior, f(c1::Coor, c2::Coor)
## n_sample:     number of MCMC sample
## burn_in:      number of samples to remove as burn-in
## n_approx:   number of Coor along one dimension to represent a domain

function runMCMC{T<:Signal}(signals::Vector{T},
                            prior_mean::Function, prior_cov::Function,
                            n_samples::Integer, burn_in::Integer=0, thinning=100,
                            n_approx::Integer=5)

    ## -----------
    ## 1) set-up

    ## create overloaded function of Prior
    f_mu, f_cov = overload_GP_function(prior_mean, prior_cov)

    ## setup coors and dicts
    coors, dic_delta_coor_index, dic_domain_coor_index = setupMCMC(signals, n_approx)

    ## Compute mean
    mu = Float64[f_mu(c) for c in coors]

    ## compute inverse covariance matrix
    Sigma = make_cov(coors, coors, f_cov)
    Sigma_inv = inv(Sigma)

    ## -----------
    ## 2) sampling

    f_sample(R) = log_joint(R, signals,
                            dic_delta_coor_index, dic_domain_coor_index,
                            mu, Sigma_inv,
                            n_approx)

    ## construct model for sampler
    mod = MCMC.model(f_sample, init=zeros(size(coors,1)))

    ## sample
    sampler_algo =  MCMC.RAM(0.1, 0.234) # RWM(0.1)
    chain = MCMC.run(mod, sampler_algo, MCMC.SerialMC(steps=n_samples, burnin=burn_in, thinning=100))

    ## -----------
    ## 3) construct dictionary compatiable to sample_predictions()

    ## create empty dictionary
    Dic_sample = Dict{Location,Vector{Float64}}()

    for i in 1:size(coors, 1)
        Dic_sample[coors[i]] = chain.samples[:,i]
    end

    return(Dic_sample)
end
