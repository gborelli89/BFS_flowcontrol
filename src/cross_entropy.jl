# -----------------------------------------------------------------------
# Cross entropy seach algorithm
# -----------------------------------------------------------------------
# Required packages
using Statistics, DataFrames
include("random_init.jl")

"""
    apply_expansion(v,γ)

## Description
Function to apply bound expansion to the domain
## Arguments
- v: value
- γ: expansion factor 
"""
#apply_expansion(v,γ) = Int(round(v*γ))

"""
    find_newdomain(data::DataFrame, nbest::Int)

## Description
Find domain range for the next iteration (CE method)
## Arguments
- data: dataframe with coordinates x and y and the values measured (works only for OneCylType by now)
- nbest: the number of best values that should be taken
"""
function find_newdomain(data::DataFrame, nbest::Int) 

    sort!(data, order(:pm, rev=true))
    best = data[1:nbest,:]
    ncoord = ncol(data) - 1

    range_coords = [minimum(best[!,i]):maximum(best[!,i]) for i in 1:ncoord]

    return range_coords
end

"""
    crossentropy_criteria(new_mean::AbstractVector, new_std::AbstractVector, previous_mean::AbstractVector; mean_tol::AbstractVector, std_tol::AbstractVector)

## Description
Criteria for cross entropy search method
## Arguments
- new_mean: array with current mean coordinate values
- new_std: array with current standard deviation coordinate values
- previous_means: array with previous mean coordinate values
- mean_tol: tolerance criteria for mean variation (for each direction)
- std_tol: tolerance criteria for standard deviation (for each direction)
"""
function crossentropy_criteria(new_mean::AbstractVector, new_std::AbstractVector, previous_mean::AbstractVector; 
                                mean_tol::AbstractVector, std_tol::AbstractVector)

    Δmean = abs.(new_mean - previous_mean)
    if all(Δmean .< mean_tol) && all(new_std .< std_tol)
        return true
    else
        return false
    end
end


"""
    crossentropy(df::AbstractVector, nbest::Int; mean_tol::AbstractVector, std_tol::AbstractVector)

## Description
Implementation of one iteration for the cross entropy seach method
## Arguments
- df: vector with all the dataframes obtained during the process (check experiments.jl)
- nbest: number of values which should be selected for the CE
- mean_tol: tolerance criteria for mean variation (for each direction)
- std_tol: tolerance criteria for standard deviation (for each direction)
## Output
Returns the bounds for a new domain.
"""
function crossentropy(df::AbstractVector, nbest::Int; mean_tol::AbstractVector, std_tol::AbstractVector) 

    n = length(df)
    i = n - 1

    if i == 0
        converge = false
    else
        ncoord = ncol(df[end]) - 1
        previous_mean = [mean(df[i][:,q]) for q in 1:ncoord]
        new_mean = [mean(df[end][:,q]) for q in 1:ncoord]
        new_std = [std(df[end][:,q]) for q in 1:ncoord]
        converge = crossentropy_criteria(new_mean, new_std, previous_mean, mean_tol=mean_tol, std_tol=std_tol)
    end

    lims = find_newdomain(df[end], nbest)
    
    return lims, converge
end


