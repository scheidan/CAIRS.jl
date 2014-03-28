## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: functions to sample predictions
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

## ---------------------------------
## sample from density p(R1*, R2*, ..., Rn* | R1, R2, ..., Rn, I1, ..., Ik)
##
## loc_pred:    vector with locations (Coor, Domain) for predictions
## R_dict_cal:  dictionary containg the samples of calibration
## n_samples:   number of samples
## prior_mean:     mean function of Prior, f(c::Coor)
## prior_cov:      covariance function of Prior, f(c1::Coor, c2::Coor)
## optional arguments:
## block_size:  number of sample points per block.
##
## Note:
## If 'block_size' is small as the number of prediction
## locations, only marginals are sampled, i.e. not all dependecies are
## visible in the sample!  Set to 'Inf' to sample a realization of the
## rain field (might be much slower).
##
## Scaling: O(n) with number of prediction points
##          O(n^3) with n_calib and blocksize

function sample_preditions{T<:Location}(loc_pred::Vector{T},
                                        R_dict_cal::Dict{Location, Vector{Float64}},
                                        n_samples::Integer,
                                        prior_mean::Function, prior_cov::Function;
                                        block_size::Real = 200)

    ## separate locations that have already been used used for calibration
    loc_pred_cal = filter(x -> in(x, collect(keys(R_dict_cal))),  loc_pred)
    loc_pred = filter(x -> !in(x, collect(keys(R_dict_cal))),  loc_pred)

    ## number of blocks
    block_size == Inf ? block_size = size(loc_pred,1) : nothing # only one block
    n_blocks = iceil(size(loc_pred,1)/block_size)

    ## initialize empty dict
    R_dict_pred = Dict{Location,Vector{Float64}}()

    ## sample each block
    for i in 1:n_blocks
        index = ((i-1)*block_size+1) : min(i*block_size, size(loc_pred, 1))
        merge!(R_dict_pred,
               sample_preditions_block(loc_pred[index], R_dict_cal, n_samples, prior_mean, prior_cov))
    end

    ## -----------
    ## predictions for non-calibration points

    n_calib = length(collect(values(R_dict_cal))[1]) # number of sample in R_dict_cal

    if length(loc_pred_cal)>0
        for loc in loc_pred_cal
            rand_index = rand(1:n_calib, n_samples)
            R_dict_pred[loc] = R_dict_cal[loc][rand_index]
        end
    end

    return(R_dict_pred)
end



## ---------------------------------
## sample from density log(p(R1*, R2*, ..., Rn* | R1, R2, ..., Rn))
##
## samples all loc_pred at ones (choleski decomposition)
##
## loc_pred:  vector with locations for predictions
## R_dict_cal:  dictionary containg the samples of calibration
## n_samples:   number of samples
## prior_mean:     mean function of Prior, f(c::Coor)
## prior_cov:      covariance function of Prior, f(c1::Coor, c2::Coor)
##
## returns a dictionary

function sample_preditions_block{T<:Location}(loc_pred::Vector{T},
                                              R_dict_cal::Dict{Location, Vector{Float64}},
                                              n_samples::Integer,
                                              prior_mean::Function, prior_cov::Function)

    loc_c = collect(keys(R_dict_cal))       # locations of calib

    ## create overlaoded function of Prior
    f_mu, f_cov = overload_GP_function(prior_mean, prior_cov)

    ## compute means
    mu_c = Float64[f_mu(loc) for loc in loc_c]
    mu_p = Float64[f_mu(loc) for loc in loc_pred]

    ## compute conditional covariance matrix
    Sigma_cc_inv = inv(make_cov(loc_c, loc_c, f_cov))
    Sigma_pp = make_cov(loc_pred, loc_pred, f_cov)
    Sigma_pc = make_cov(loc_pred, loc_c, f_cov)

    Sigma_temp = Sigma_pc * Sigma_cc_inv
    Sigma_cond = Sigma_pp - Sigma_temp * Sigma_pc'
    Sigma_cond_chol = chol(Sigma_cond, :L)

    ## number of sample in R_dict_cal
    n_calib = length(collect(values(R_dict_cal))[1])
    n_pred = size(loc_pred, 1)

    ## Array to store the samples temporary
    R_array_pred = Array(Float64, n_samples, n_pred)

    for i in 1:n_samples

        ## pick a random realization
        rand_index = rand(1:n_calib)
        R = Float64[]
        for c in keys(R_dict_cal)
            push!(R, R_dict_cal[c][rand_index])
        end

        ## mean
        mean = mu_p + Sigma_temp * (R - mu_c)

        ## sample from multivariate normal
        R_array_pred[i,:] = mean + Sigma_cond_chol*randn(n_pred)

    end

    ## write 'R_array_pred' in Dictionary
    R_dict_pred = Dict{Location,Vector{Float64}}()
    for i in 1:length(loc_pred)
        R_dict_pred[loc_pred[i]] = R_array_pred[:,i]
    end

    return(R_dict_pred)
end
