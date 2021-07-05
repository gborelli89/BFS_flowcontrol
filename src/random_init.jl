# -----------------------------------------------------------------------
# Initiate population randomly accordingly to a distribution
# -----------------------------------------------------------------------
# Required packages
using Interpolations, StatsBase
include("stepper_move.jl")

"""
    randweights(coord::AbstractArray, p::AbstractVector, α = 1.0)

## Description
Find the weights for a random sampling considering one dimension only (used for the upstream cylinder position)
## Arguments
- coord: coordinate array
- p: pressure array measured at coord
- α: magnification factor (p^α)
## Output
A weight function (linearly interpolated) is returned
## Example
```jldoctest
julia> x = 0:15:90;
julia> p = [-95,-90,-80,-72,-80,-90,-95];

julia> itpx = randweights(x, p, α=2.0);
julia> itpx(50)
0.4156138645934563
```
"""
function randweights(coord::AbstractArray, p::AbstractVector; α = 1.0)

    val = p .- maximum(p)
    val = val .- minimum(val)
    prob = val.^α ./ sum(val.^α)
    itp = LinearInterpolation(coord, prob)

    return itp
end


"""
    randpos(n::Int, lims::AbstractArray)

## Description
Function to generate random position of OneCylType
## Arguments
- n: number of positions
- lims: array with limits   
- weights: array with probabilities for the upstream cylinder position (x and y positions) 
## Example
```jldoctest
julia> p = randompos(5, [0:90,-25:25])
5-element Array{OneCylType,1}:
 OneCylType(50.0, -23.0)
 OneCylType(76.0, -22.0)
 OneCylType(85.0, 18.0)
 OneCylType(48.0, 15.0)
 OneCylType(6.0, 25.0)
```
"""
function randpos(n::Int, lims::AbstractArray)

    if length(lims) == 2
        pos = OneCylType.(sample(lims[1],n), sample(lims[2],n))
    elseif length(lims) == 4
        pos = TwoCylType.(sample(lims[1],n), sample(lims[2],n), sample(lims[3],n), sample(lims[4],n))
    elseif length(lims) == 6
        pos = ThreeCylType.(sample(lims[1],n), sample(lims[2],n), sample(lims[3],n), sample(lims[4],n),
                            sample(lims[5],n), sample(lims[6],n))
    else
        throw(BoundsError)
    end
end

function randpos(n::Int, lims::AbstractArray, weigths::AbstractArray)

    if length(lims) == 2
        pos = OneCylType.(sample(lims[1],Weights(weigths[1]),n), sample(lims[2],Weights(weigths[2]),n))
    elseif length(lims) == 4
        pos = TwoCylType.(sample(lims[1],Weights(weigths[1]),n), sample(lims[2],Weights(weigths[2]),n), 
                        sample(lims[3],n), sample(lims[4],n))
    elseif length(lims) == 6
        pos = ThreeCylType.(sample(lims[1],Weights(weigths[1]),n), sample(lims[2],Weights(weigths[2]),n), 
                        sample(lims[3],n), sample(lims[4],n), sample(lims[5],n), sample(lims[6],n))
    else
        throw(BoundsError)
    end
end

"""
    randindividual(nc=3; size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4], d=5.0)

## Description
Random initialization
## Arguments
- nc: number of cylinders
- size: size of the search 2^size (mm)
- coord_zero: initial coordinates (offset)
- d: cylinder diameter + offset (mm)
## Example
```jldoctest
julia> randindividual(1)
```
"""
function randindividual(nc=3; size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4], d=5.0)
    xs = []
    ys = []
    
    if length(size) != length(coord_zero)
        error("size and coord_zero must have the same length!")
    end

    for i in 1:2:length(size)
        push!(xs, rand(0:2^(size[i])-1) + coord_zero[i])
        k = i+1
        push!(ys, rand(0:2^(size[k])-1) + coord_zero[k])
    end

    if nc == 1 
        p = OneCylType(xs[1], ys[1])
    elseif nc == 2
        p = TwoCylType(xs[1], ys[1], xs[2]+d, ys[2])
    elseif nc == 3
        p = ThreeCylType(xs[1], ys[1], xs[2]+d, ys[2], xs[3]+d, ys[3])
    else
        error("Number of cylinders not defined")
    end
    return p
end

"""
    population_init(n, nc=3; size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4], d=5.0)

## Description
Create initial population with random individuals.
## Arguments
- n: number os individuals
- nc: number of cylinders
- size: size of the search 2^size (mm)
- coord_zero: initial coordinates (offset)
- d: cylinder diameter + offset (mm)
## Example
```jldoctest
# Population with 2 control cylinders and 20 individuals
julia> p = population_init(20,2);
```
"""
function population_init(n, nc=3; size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4], d=5.0)
    P = [randindividual(nc, size=size, coord_zero=coord_zero, d=d) for _ in 1:n]
end


