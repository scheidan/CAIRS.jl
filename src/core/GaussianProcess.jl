## =======================================================
## Continuous Assimilation of Integrating Rain Sensors (CAIRS)
##
## Description: Gaussian process (GP) as prior for rain field
##
## Andreas Scheidegger -- andreas.scheidegger@eawag.ch
## =======================================================


## ---------------------------------
## transformation
##
## if X~N(0,1), g(X) should have a realistic rain distibution

## exponential with offset
## function trans2real(R_trans::Float64)
##     k = 0.8
##     offset = 1.5
##     R_trans <= offset ? 0 : exp(k*(R_trans-offset)) - 1.0
## end

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

## ## vectorized version
## function trans2real(R_trans::Vector{Float64})
##     n = size(R_trans,1)
##     R_real = zeros(n)
##     for i in 1:n
##         R_real[i] = trans2real(R_trans[i])
##     end
##     return(R_real)
## end

## ---------------------------------
## function that construct overloaded mean and covariance functions

function overload_GP_function(f_mean::Function, f_covariance::Function)

    !method_exists(f_mean, (Coor,)) ?
    error("The mean function of the GP must provide a method for arguments of type 'Coor'!") : nothing

    !method_exists(f_covariance, (Coor,Coor)) ?
    error("The covariance function of the GP must provide a method for both arguments of type 'Coor'!") : nothing

    ## ---------------------------------
    ## mean function

    f_mu(c::Coor) = f_mean(c)

    ## --- for Domains

    function f_mu(d::Domain)

        ## vector of domain extend
        extend_v = [d.extend.x,
                    d.extend.y,
                    d.extend.time]

        ## vector of positions
        pos_v = [d.position.x,
                 d.position.y,
                 d.position.time]

        ## construct function to integrate over
        function f_int(v::Vector{Float64})
            ## change the coordinates for integration
            kk = 1
            for i in 1:3
                if extend_v[i] != 0.0
                    pos_v[i] = v[kk]
                    kk += 1
                end
            end

            f_mu(rotate(Coor(pos_v[1], pos_v[2], pos_v[3]), d.position, d.angle))
        end

        ## construct lower and upper bounds
        lower = pos_v[extend_v.!=0.0]
        upper = (pos_v + extend_v)[extend_v.!=0.0]

        ## integrate
        I = pcubature(f_int, lower, upper, maxevals=10000)
        ## println(I, ", relative error: ", I[2]/I[1])
        return(I[1])
    end


    ## ---------------------------------
    ## covariance function

    f_cov(c1::Coor, c2::Coor) = f_covariance(c1, c2)

    ## -----------
    ## covariance (kernel) function for ( Domain, Domain )

    function f_cov(d1::Domain, d2::Domain)

        ## vector of domain extend
        extend_v = [d1.extend.x,
                    d1.extend.y,
                    d1.extend.time,
                    d2.extend.x,
                    d2.extend.y,
                    d2.extend.time]

        ## vector of positions
        pos_v = [d1.position.x,
                 d1.position.y,
                 d1.position.time,
                 d2.position.x,
                 d2.position.y,
                 d2.position.time]


        ## construct function to integrate over
        function f_int(v::Vector{Float64})
            ## change the coordiates for integration
            kk = 1
            for i in 1:6
                if extend_v[i] != 0.0
                    pos_v[i] = v[kk]
                    kk += 1
                end
            end

            f_cov(rotate(Coor(pos_v[1], pos_v[2], pos_v[3]), d1.position, d1.angle),
                  rotate(Coor(pos_v[4], pos_v[5], pos_v[6]), d2.position, d2.angle))
        end

        ## construct lower and upper bounds
        lower = pos_v[extend_v.!=0.0]
        upper = (pos_v + extend_v)[extend_v.!=0.0]

        ## integrate
        if sum(extend_v .!= 0.0)>2          # for >2-dimensional integration
            I = hcubature(f_int, lower, upper, maxevals=10000)
        else
            I = pcubature(f_int, lower, upper, maxevals=10000)
        end
        ## println(I, ", relative error: ", I[2]/I[1])
        return(I[1])

    end

    ## --- for ( Coor, Domain )

    function f_cov(c::Coor, d::Domain)
        d1 = Domain(c, Coor(0,0,0), 0.0)    # a coordinate is a domain with a zero extend
        f_cov(d1, d)
    end

    f_cov(d::Domain, c::Coor) = f_cov(c, d)

    return(f_mu, f_cov)                 # retun overloaded functions


end

## ---------------------------------
## construct covariance matrix


function make_cov{T1<:Location, T2<:Location}(loc_1::Vector{T1}, loc_2::Vector{T2},
                  f_cov::Function)
    ## if asymmetric
    if loc_1 != loc_2
        Sigma = Float64[f_cov(l1, l2) for l1 in loc_1, l2 in loc_2]
    else
        ## for symmetric matrices
        n = length(loc_1)
        Sigma = zeros(n, n)
        @inbounds for j in 1:n
            for i in j:n
                if i !=j
                    Sigma[i,j] = Sigma[j,i] = f_cov(loc_1[i], loc_2[j])
                else
                    Sigma[i,j] = f_cov(loc_1[i], loc_2[j])
                end
            end
        end

    end

    return(Sigma)
end




## ---------------------------------
## proportional to joint density p(R1, R2, ..., Rn, I1, I2, ...)

function log_p_prior{T<:Location}(locations::Vector{T}, Samp_dict::Dict{Location, Vector{Float64}},
                     i_sample::Int,
                     mu::Vector{Float64},
                     Sigma::PDMats.AbstractPDMat)

    ## get rain at all locations
    R = Float64[]
    for loc in locations
        push!(R, Samp_dict[loc][i_sample])
    end

    ## log density
    D = R - mu
    neg_log_p = PDMats.invquad(Sigma, D) # compute D' inv(Sigma) * D
    return(-neg_log_p)
end
