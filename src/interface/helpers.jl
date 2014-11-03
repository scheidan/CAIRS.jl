## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: little helper functions
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================


## ---------------------------------
## Write MCMC chains as csv file

## samples_dict: dictionary holding the MCMC sampler of the calibration points
## filename:     string of filename

function chains2csv(samples_dict::Dict{Location,Vector{Float64}}, filename="chains.csv")

    loc1 = collect(keys(samples_dict))[1]
    n_sample = size(samples_dict[loc1], 1)

    ## write Samp_dict as array
    chains = Array(Float64, n_sample, length(samples_dict))
    for i in 1:length(samples_dict)
        loc = collect(keys(samples_dict))[i]
        chains[:,i] = samples_dict[loc]
    end
    writecsv(filename, chains)
end

## ---------------------------------
## Write summary of predicted loactions

## pred_dict:   dictionary of predicted Coor's
## filename:     string of filename

function summary2csv(pred_dict::Dict, filename="predictions.csv", real::Bool=true)

    n_coor = mapreduce(x -> typeof(x)==Coor, +, keys(pred_dict))
    predictions = Array(Float64, n_coor, 7)

    for i in 1:n_coor
        coor = collect(keys(pred_dict))[i]
        if typeof(coor)==Coor
            Rcoor = pred_dict[coor]
            predictions[i,:] = [coor.x, coor.y, coor.time,
                          mean(Rcoor),
                          std(Rcoor),
                          quantile(Rcoor, 0.1),
                          quantile(Rcoor, 0.9)]'
        end
    end

    writecsv(filename, predictions)
end



## ---------------------------------
## Write coordinates and domains of sensors as csv

## signal_dict: vector of signals
## filename:    string of filename

function sensor2csv(signals::Vector, filename="sensors.csv")

    positions = Array(Any, 0, 4)

    x = Float64[]
    y = Float64[]
    time = Float64[]
    name = ASCIIString[]

    d = 1
    for sig in signals

        ## Coordinates
        if length(sig.sensor.delta_coor) > 0
            for c in sig.sensor.delta_coor
                coor = sig.position + c
                push!(x, coor.x)
                push!(y, coor.y)
                push!(time, coor.time)
                push!(name, "coor")
            end
        end

        ## integration domains, write coordinates of each corner
        if sig.sensor.domain_extent != Coor(0.0, 0.0, 0.0)

            for b1 in 0:1
                for b2 in 0:1
                    for b3 in 0:1
                        ## shift and rotate coordiante
                        coor = rotate(Coor(sig.position.x + sig.sensor.domain_extent.x*b1,
                                           sig.position.y + sig.sensor.domain_extent.y*b2,
                                           sig.position.time + sig.sensor.domain_extent.time*b3),
                                      sig.position, sig.angle)

                        push!(x, coor.x)
                        push!(y, coor.y)
                        push!(time, coor.time)
                        push!(name, "domain_$d")
                    end
                end
            end
            d += 1
        end
    end

    ## write file
    writecsv(filename, [x y time name])
end
