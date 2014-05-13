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
## n_approx:   number of points per dimension to approximate domains

function setupMCMC{T<:Signal}(signals::Vector{T}, n_approx::Integer=5)

    coors = Coor[]                      # array holds all coordinates that are 'measured'
    dic_delta_coor_index = [S => Integer[] for S in signals]  # Signal => [indices of coors directly measured]
    dic_domain_coor_index = [S => Integer[] for S in signals]  # Signal => [indices of coors for approx. domains]

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
        ## approximate domains with a points

        if S.sensor.domain_extent != Coor(0.0, 0.0, 0.0)
            ## rotate and shift in right position
            coors_approx = unique([rotate(Coor(x, y, time), Coor(0.0, 0.0, 0.0), S.angle) + S.position
                                   for x in linspace(0, S.sensor.domain_extent.x, n_approx),
                                       y in linspace(0, S.sensor.domain_extent.y, n_approx),
                                       time in linspace(0, S.sensor.domain_extent.time, n_approx)])
            coors_approx = reshape(coors_approx, length(coors_approx)) # make vector


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
## computes log likelihood of a signals

## R:        vector with rain intensities corresponding to coors
## signals:  vector of 'Signal' object
## dic_delta_coor_index:     dictionaries as created by setupMCMC()
## dic_domain_coor_index:    dictionaries as created by setupMCMC()

function log_likeli{T<:Signal}(R::Vector{Float64}, signals::Vector{T},
                               dic_delta_coor_index, dic_domain_coor_index)

    ll = 0
    for S in signals

        if S.sensor.domain_extent == Coor(0.0, 0.0, 0.0) # no domain
            ll += S.sensor.log_p(S.signal, R[dic_delta_coor_index[S]])
        end
        if size(S.sensor.delta_coor,1) ==  0 # no coordiantes
            ll += S.sensor.log_p(S.signal, mean(R[dic_domain_coor_index[S]]))
        end
        if (S.sensor.domain_extent != Coor(0.0, 0.0, 0.0)) && (size(S.sensor.delta_coor,1) > 0)
            ll += S.sensor.log_p(S.signal, R[dic_delta_coor_index[S]], mean(R[dic_domain_coor_index[S]]))
        end

    end
    return(ll)
end
