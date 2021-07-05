# -----------------------------------------------------------------------
# Experiments pre processing
# -----------------------------------------------------------------------
# Required packages
using JLD, DataFrames, Dates, Statistics

"""
    function op2column(A::AbstractArray, v::AbstractVector, op=+)

## Description
Apply operation for the columns of a matrix
## Arguments
- A: AbstractArray
- v: vector to perform elementwise operation
- op: operation
## Example
```jldoctest
julia> A = [[1,2,3] [4,5,6] [7,8,9]]
3×3 Array{Int64,2}:
 1  4  7
 2  5  8
 3  6  9

julia> v = [1,1,1]
 3-element Array{Int64,1}:
 1
 1
 1

julia> op2column(A, v, op=-)
3×3 Array{Int64,2}:
 0  3  6
 1  4  7
 2  5  8
```
"""
function op2column(A::AbstractArray, v::AbstractVector; op=+)
    return op.(A, reshape(repeat(v, inner=size(A)[1]), size(A)))
end

"""
    correct_zero(fname::String, fzero_init::String, fzero_end::String)

## Decription
Correct zero by interpolating zero measurement at the begining and at the end.
## Arguments
- fname: string (or array) with the name of the test file
- fzero_init: name of the file with the zero read before the test
- fzero_end: name of the file with the zero read after the test 
"""
function correct_zero(fname::String, fzero_init::String, fzero_end::String)

    p,h = load(fname, "pressure", "h")
    pz1,hz1 = load(fzero_init, "pressure", "h")
    pz2,hz2 = load(fzero_end, "pressure", "h") 

    Δhz = (hz2 - hz1)/(Millisecond(1.0)*1000) # dt in seconds
    Δh = (h .- hz1) ./ (Millisecond(1.0)*1000)

    pz1_m = [mean(pz1[:,i]) for i in 1:16]
    pz2_m = [mean(pz2[:,i]) for i in 1:16]
    pz_m = pz2_m .- pz1_m

    Δp = [pz1_m + (x/Δhz)*pz_m for x in Δh]

    p_corr = [op2column(p[i], Δp[i], op=-) for i in 1:length(p)]

    return p_corr
end

function correct_zero(fname::AbstractArray, fzero_init::String, fzero_end::String)

    nf = length(fname)
    pz1,hz1 = load(fzero_init, "pressure", "h")
    pz2,hz2 = load(fzero_end, "pressure", "h") 
    
    pz1_m = [mean(pz1[:,i]) for i in 1:16]
    pz2_m = [mean(pz2[:,i]) for i in 1:16]
    pz_m = pz2_m .- pz1_m
    
    Δhz = (hz2 - hz1)/(Millisecond(1.0)*1000) # dt in seconds

    p = []
    h = []
    for i in 1:nf
        ptemp,htemp = load(fname[i], "pressure", "h")
        push!(p, ptemp)
        push!(h, htemp)
    end
    
    Δh = (h .- hz1) ./ (Millisecond(1.0)*1000)
    Δp = [pz1_m + (x/Δhz)*pz_m for x in Δh]

    p_corr = [op2column(p[i], Δp[i], op=-) for i in 1:length(p)]

    return p_corr
end


"""
    calc_meanpressure(p::AbstractArray, nch)

## Description
Compute the mean pressure values given an array of scanivalve measurements
## Arguments
- p: pressure measured data (already corrected if desired)
- nch: scanivalve channels to consider
"""
function calc_meanpressure(p::AbstractArray, nch)
    pm = [mean(i[:,nch]) for i in p]
end

"""
    calc_meancp(p::AbstractArray, nch)

## Description
Compute the mean pressure coefficient values given an array of scanivalve measurements
## Arguments
- p: pressure measured data (already corrected if desired)
- nch: scanivalve channels to consider
- chref: reference channel (dinamic pressure measurement)
"""
function calc_meancp(p::AbstractArray, nch, chref=16)
    pref = [mean(i[:,chref]) for i in p]
    cpm = calc_meanpressure(p, nch) ./ pref
    return [cpm, pref]
end


"""
    calc_stdpressure(p::AbstractArray, nch)

## Description
Compute the standard deviation pressure values given an array of scanivalve measurements
## Arguments
- p: pressure measured data (already corrected if desired)
- nch: scanivalve channels to consider
"""
function calc_stdpressure(p::AbstractArray, nch)
    pstd = [std(i[:,k]) for i in p for k in nch]
    pstd = reshape(pstd, (length(nch), length(p)))
    return pstd
end

"""
    calc_stdcp(p::AbstractArray, nch)

## Description
Compute the standard deviation pressure coefficient values given an array of scanivalve measurements
## Arguments
- p: pressure measured data (already corrected if desired)
- nch: scanivalve channels to consider
- chref: reference channel (dinamic pressure measurements)
"""
function calc_stdcp(p::AbstractArray, nch, chref=16)
    cpstd = [std(i[:,k] ./ i[:,chref]) for i in p for k in nch]
    cpstd = reshape(cpstd, (length(nch), length(p)))
    return cpstd
end

