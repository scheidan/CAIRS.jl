## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: calibration of rainfield
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

## ---------------------------------
## Generates a dictionary of all signals
##
## Assign all signals that are influenced
## by the intensity to a coordinate or a domain:
## Domain/Coor => [signal1, Signal3, ...]

## signals:    vector of Signals

function make_signal_dict{T<:Signal}(signals::Vector{T})

    ## create empty dictionary
    Dic = Dict{Location, Vector{Signal}}()

    ## loop over all Signals
    for S in signals

        ## construct integration domain and add to Dict
        if S.sensor.domain_extent != Coor(0.0, 0.0, 0.0)
            domain = Domain(S.position, S.sensor.domain_extent, S.angle)
            ## if 'domain' is already a key in Dic
            if haskey(Dic, domain)
                Dic[domain] = unique([S, Dic[domain]])
                ## if 'domain' is a new key
            else
                Dic[domain] = [S]
            end
        end

        ## construct all coordinated on that S is conditioned
        if size(S.sensor.delta_coor, 1) > 0
            coors = Coor[]
            for delta_coor in S.sensor.delta_coor
                push!(coors, S.position + delta_coor)
            end

            ## add all coors to Dict
            for i in 1:size(coors, 1)
                ## if coors[i] is already a key in Dic
                if haskey(Dic, coors[i])
                    Dic[coors[i]] = unique([S, Dic[coors[i]]])
                    ## if coors[i] is a new key
                else
                    Dic[coors[i]] = [S]
                end
            end
        end

    end
    return(Dic)

end

## ---------------------------------
## Generates a dictionary of all sample points
##
## Assign to each coordinate 'Coor' or 'Domain" a vector for
## for MCMC samples.
## Coor/Domain => [R_1, R_2, ..., R_{n_samples}]

## signal_dict:   dictionary of signals
## n_sample:      length of the initialized vector for MCMC sample

function make_sample_point_dict(signal_dict::Dict{Location, Vector{Signal}},
                                n_samples::Int)

    ## create empty dictionary
    Dic = Dict{Location,Vector{Float64}}()

    ## use the same keys as for signal.dict
    for c in keys(signal_dict)
        Dic[c] = zeros(n_samples)
    end
    return(Dic)
end


## ---------------------------------
## computes log_p of a signal given a sample point dictionary

## S:           a 'Signal' object
## sample_dict: dictionary of all sample points
## i_sample:    index of MCMC sample

function log_p_of_signal(S::Signal, sample_dict::Dict{Location, Vector{Float64}}, i_sample::Int)

    ## find all coordinates on that S is conditioned
    if size(S.sensor.delta_coor,1) > 0
        coors = Coor[]
        for delta_coor in S.sensor.delta_coor
            push!(coors, S.position + delta_coor)
        end

        ## get vector of all intensities
        R = Float64[]
        for c in coors
            push!(R, sample_dict[c][i_sample])
        end
    end


    ## get value of integrated Domain
    if S.sensor.domain_extent != Coor(0.0, 0.0, 0.0)
        domain = Domain(S.position, S.sensor.domain_extent, S.angle)
        I = sample_dict[domain][i_sample]
    end

    ## compute log_p
    if S.sensor.domain_extent == Coor(0.0, 0.0, 0.0) # no domain
        log_p = S.sensor.log_p(S.signal, R)
    end
    if size(S.sensor.delta_coor,1) ==  0 # no coordiantes
        log_p = S.sensor.log_p(S.signal, I)
    end
    if (S.sensor.domain_extent != Coor(0.0, 0.0, 0.0)) && (size(S.sensor.delta_coor,1) > 0)
        log_p = S.sensor.log_p(S.signal, R, I)
    end
    return(log_p)
end


## ---------------------------------
## Gibbs sampling with adaptive jump distributions
##
## Roberts G, Rosenthal J (2009). "Examples of Adaptive MCMC."
##   Computational Statistics and Data Analysis, 18, 349-367.
##
## signals:      vector of 'Signals'
## prior_mean:     mean function of Prior, f(c::Coor)
## prior_cov:      covariance function of Prior, f(c1::Coor, c2::Coor)
## n_sample:     number of MCMC sample
## burn_in:      number of samples to remove as burn-in
## adaption:     true/false

## FUTURE: - Ordered Overrelaxation (see MacKay)?
##         - adaptive rejection sampling (http://www.stat.duke.edu/~cnk/Links/slides.pdf)?
##         - slice sampling?

function Gibbs{T<:Signal}(signals::Vector{T},
                          prior_mean::Function, prior_cov::Function,
                          n_samples::Int, burn_in::Int=0;
                          adaption::Bool=true)

    ## -----------
    ## 1) set-up

    ## create overlaoded function of Prior
    f_mu, f_cov = overload_GP_function(prior_mean, prior_cov)

    ##  create dictionary of all signals
    signal_dict = make_signal_dict(signals)

    ## create dictionary of all point rains to sample from
    samples_dict = make_sample_point_dict(signal_dict, n_samples)
    samples_dict_prop = make_sample_point_dict(signal_dict, 1)

    ## array with all coordiates
    samp_points = collect(keys(samples_dict))

    ## Compute mean
    mu = Float64[f_mu(loc) for loc in samp_points]

    ## compute inverse covariance matrix
    Sigma = PDMats.PDMat(make_cov(samp_points, samp_points, f_cov))

    ## sd of jump distributions
    sd_prior = sqrt.(diag(Sigma))
    sd_prop = Dict(collect(keys(samples_dict))[ii] => sd_prior[ii] for ii in 1:length(samples_dict))

    ## -----------
    ## 2) sampling

    for i in 1:n_samples-1

        mod(i, 5000)==0 ? println("iteration $i of $n_samples") : nothing

        ## loop over all all point rains
        for location in keys(samples_dict)

            ## find all signals that condition on 'location'
            signals = signal_dict[location]

            ## proposal rain
            samples_dict_prop[location][1] = samples_dict[location][i] + randn()*sd_prop[location]

            ## compute acceptance probability
            log_p_prop = log_p_prior(samp_points, samples_dict_prop, 1, mu, Sigma)
            for S in signals
                log_p_prop += log_p_of_signal(S, samples_dict_prop, 1)
            end

            log_p_old = log_p_prior(samp_points, samples_dict, i, mu, Sigma)
            for S in signals
                log_p_old += log_p_of_signal(S, samples_dict, i)
            end

            p_acc = min(1.0, exp(log_p_prop - log_p_old))

            ## accept or not
            if rand() < p_acc[1]
                samples_dict[location][i+1] = samples_dict_prop[location][1]
            else
                samples_dict[location][i+1] =  samples_dict[location][i]
                samples_dict_prop[location][1] = samples_dict[location][i]
            end


            ## adapt jump distribution, see Roberts G, Rosenthal J (2009), Sec. 3
            ## (adaption is only during burn-in active)
            if adaption && i < burn_in && mod(i, 50)==0
                acc_rate = length(unique(samples_dict[location][1:i+1]))/(i+1) # current acceptance rate
                if acc_rate > 0.44
                    sd_prop[location] = sd_prop[location] * exp(min(0.1, 1/sqrt(i)))
                else
                    sd_prop[location] = sd_prop[location] / exp(min(0.1, 1/sqrt(i)))
                end

            end

        end

    end

    ## -----------
    ## 3) remove burn-in samples

    if burn_in > 0
        for location in keys(samples_dict)
            splice!(samples_dict[location], 1:burn_in)
        end
    end

    return(samples_dict)
end
