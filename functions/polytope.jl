ENV["R_HOME"] = "C:\\Program Files\\R\\R-4.3.2"

using RCall 

# install Rcpp in R, open R terminal and type: install.packages("Rcpp"), do the same for volesti: install.packages("volesti")
R"library(Rcpp)" # uses Rcall package to run R code in julia -> make sure R is installed on your device!


function create_pol(N::Integer)
    """Create N-dim unit cube."""
    rSTDO = R"""
    library(volesti)
    P <- gen_cube($N,'H')
    values_b <- as.numeric(cbind(matrix(0,1,$N), matrix(1,1,$N))) 
    slot(P,'b') <- values_b
    """ # mondify cube from [-1 1] to [0 1]
    return rcopy(R"slot(P, 'A')"),rcopy(R"slot(P, 'b')")
end

function comp_vol()
    R"""
    v <- volume(P, settings = list("error" = 10^-12))
    """
    return rcopy(R"v")
end

function return_P()
    return rcopy(R"P")
end

#R"""
#library(volesti)
#P = gen_rand_hpoly(20, 60, generator = list('constants' = 'sphere'))"""

function create_scaled_pol(N::Integer,L::AbstractArray,U::AbstractArray)
    """Create polytope with increased upper bounds."""
        rSTDO = R"""
        library(volesti)
        P <- gen_cube($N,'H')
        values_b <- as.numeric(cbind(matrix(-c(unlist($L)),1,$N), matrix(c(unlist($U)),1,$N)))
        slot(P,'b') <- values_b
        #Hpolytope$new(P$A, P$b)
        """ # mondify cube from [-1 1] to [0 1]
        return rcopy(R"slot(P, 'A')"),rcopy(R"slot(P, 'b')"),rcopy(R"P")
end


function sample_pol(number_of_samples::Integer=1)
    samples = rcopy(R"sample_points(P, $number_of_samples)")
    return samples::Array{Float64,2}
end

function sample_pol_Walk(number_of_samples::Integer=1; random_walk::String="BRDHR")
    rSTDO = R"""
    library(volesti)
    # Create a list for the random walk parameter
    random_walk_list <- list("walk" = $random_walk)
    # Sample points
    samples <- sample_points(P, n = $number_of_samples, list("walk" = $random_walk))
    """
    samples = rcopy(R"samples")
    return samples
end

function get_pol()
    A = rcopy(R"slot(P, 'A')")
    b = rcopy(R"slot(P, 'b')")
    return A::Array{Float64,2}, b::Union{AbstractArray{<:Number},Number}
end


function add_to_pol(Al::AbstractArray,bl::Union{AbstractArray{<:Number},Number})
    R"""
    slot(P, 'A') <- rbind(slot(P, 'A'), $Al)
    slot(P, 'b') <- c(slot(P, 'b'), $bl)"""
    #nothing
    return rcopy(R"P")
end


function get_full_pol()
    A = rcopy(R"slot(P, 'A')")
    b = rcopy(R"slot(P, 'b')")
    P = rcopy(R"P")
    return P, A::Array{Float64,2}, b::Union{AbstractArray{<:Number},Number}
end

function remove_from_pol()
    R"""
    # Get the current 'A' and 'b' slots
    A <- slot(P, 'A')
    b <- slot(P, 'b')
    
    # Remove the last row from 'A' and the last element from 'b'
    if (nrow(A) > 0 && length(b) > 0) {
        slot(P, 'A') <- A[-nrow(A), , drop = FALSE]
        slot(P, 'b') <- b[-length(b)]
    }
    """
    return rcopy(R"P")
end









