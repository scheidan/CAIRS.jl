## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: high-level interface for prediction
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================

## ---------------------------------
## Predict (integrated) intensities condition on all relevant signals
##
## loc_pred:       Vector or array of 'Location's for predictions
## signals:        Vector or array of 'Signal's
## - optional: -
## n_sample_calib: number of MCMC samples for calibration
## burn_in:        number of calibration samples that are discarded
## n_sample_pred:  number of samples for prediction
## delta:          maximal time distance of signals to prediction location in Milli seconds


function predict{T<:Location}(loc_pred::Array{T}, signals::Vector{Signal};
                 n_sample_calib::Int = 20000, burn_in::Int = -1,
                 n_sample_pred::Int = 5000, delta::Real = -1.0)

    ## reshape to Vector
    loc = reshape(loc_pred, (length(loc_pred),))

    ## select the relevant Signals
    delta == -1.0 ? delta = 3 * l_temporal : nothing
    sig = find_near_signals(loc, signals, delta)

    if size(sig,1) == 0
        println("No signals are close enough. Increase 'delta'.")
        return()
    end

    println("\n--- calibration ---")
    println(" Consider ", size(sig, 1), " signals")
    Samp_dict_cal = Gibbs(sig, n_sample_calib, burn_in, true)
    
    ## sample predictions

    println("\n--- prediction ---")
    println(" ", length(loc), " locations")
    R_dict_pred = sample_preditions(loc, Samp_dict_cal, n_sample_pred, block_size=200)

    R_dict_pred

end