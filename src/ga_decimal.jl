# -----------------------------------------------------------------------
# Decimal coded genetic algorithm
# -----------------------------------------------------------------------
# Reference: JENKINS, W. M. A decimal-coded evolutionary algorithm for 
# constrained optimization. Computers and Structures, v. 80, 
# p. 471–480, 2002.
# -----------------------------------------------------------------------
# Required packages
using DataFrames, Statistics
include("random_init.jl")

# Algorithm implemented according to the following paper:
# "A decimal-coded evolutionary algorithm for constrained optimization", W.M. Jenkins, 2002
# The algorithm is decimal coded, mutation based. There is no crossover.

# Mutation

"""
    αmax(gen::Int; α0::Real, maxgen::Int)

## Description
Function to compute the maximum shift for the generation
## Arguments
- gen: current generation
- α0: shift for the first generation
- maxgen: total number of generations
"""
αmax(gen::Int; α0::Real, maxgen::Int) = α0*(1-gen/maxgen)^2 + 1

"""
αind(gen::Int, ind::Int; α0::Real, maxgen::Int, ucb::Int, popsize::Int, rounded=true)

## Description
Function to compute the shift for each individual acconding to the current generation
## Arguments
- gen: current generation
- ind: individual id
- α0: shift for the first generation
- maxgen: total number of generations
- ucb: upper class boundary (elitism)
- popsize: population size
- digits: decimal places
"""
function αind(gen::Int, ind::Int; α0::Real, maxgen::Int, ucb::Int, popsize::Int, digits=0)

    shift_gen = αmax(gen, α0=α0, maxgen=maxgen)

    α = (ind - ucb)*(shift_gen - 1)/(popsize - ucb) + 1
    α = round(α, digits=digits)

    return α
end

# Elitism

"""
    elite!(pop::DataFrame; lcb::Int, cols= :fitness)

## Description
Order population from the most well adapted induvidual to the least adapted. 
Also, remove individuals from the lower boundary and replace them by the most well adapted
## Arguments
- pop: population dataframe with coordinates and fitness
- lcb: lower class boundary
- cols: idicates the finess column
"""
function elite!(pop::DataFrame; lcb::Int, cols= :fitness)
    
    sort!(pop, cols)
    delete!(pop, 1:lcb)
    sort!(pop, cols, rev=true)
    append!(pop, pop[1:lcb,:])

    return pop
end

# Boundaries

"""
    isindomain(coord::Real, lims)

## Description
Identify if a coordinate is inside a domain
## Arguments
- coord: coordinate
- lims: boundaries (range or 1D array)
## Example
```jldoctest
julia> isindomain(10, 1:4)
false

julia> isindomain(10, 1:11)
true
```
"""
isindomain(coord::Real, lims) = coord >= lims[1] && coord <= lims[end]

"""
    findnearest(val::Real, lims)

## Description
Find the nearest neighbor in the case the coordinate is outside the domain
## Arguments
- coord: coordinate
- lims: noundaries (range or 1D array) 
"""
function findnearest(coord::Real, lims)

    if coord < lims[1]
        return lims[1]
    elseif coord > lims[end]
        return lims[end]
    else
        error("The coordinate should be outside the domain!")
    end
end

"""
    indomain_coord!(coord::Real, boundaries)

## Description
Make all the coordinates stay inside the boundary
## Arguments
- coord: coordinate
- boundaries: domain limits for the direction
"""
function indomain_coord(coord::Real, boundaries)

    if !isindomain(coord, boundaries)
        coord_new = findnearest(coord, boundaries)
    else
        coord_new = coord
    end
    
    return coord_new
end


"""
    gadecimal_step(pop::DataFrame; Tm::AbstractVector, α0::AbstractVector, gen::int, maxgen::Int, ucb::Int, lcb::Int, digits=0)

## Description
One step of decimal-coded GA algorithm
## Arguments
- pop: input population
- Tm: mutation rate for each direction
- α0: shift for the first generation
- gen: current generation
- maxgen: total number of generations
- ucb: upper class boundary (elitism)
- lcb: lower class boundary
- digits: decimal places
"""
function gadecimal_step(pop::DataFrame; Tm::AbstractVector, α0::AbstractVector, gen::Int, maxgen::Int, 
                        ucb::Int, lcb::Int, digits=0, boundaries::AbstractVector)
    
    n = nrow(pop)

    pop_new = deepcopy(pop)
    elite!(pop_new, lcb=lcb)
    new_individuals = []

    for i in (ucb+1):n
        
        ncoord = length(α0)
        
        α = [αind(gen, i, α0=a, maxgen=maxgen, ucb=ucb, popsize=n, digits=digits) for a in α0]
        mutate = rand(ncoord) .> Tm
        
        if sum(mutate) != 0
            push!(new_individuals, i)
        end

        for j in 1:ncoord
            if !mutate[j]
                nothing
            else
                op = rand([+,-])
                pop_new[i,j] = op(pop_new[i,j], α[j])
                pop_new[i,j] = indomain_coord(pop_new[i,j], boundaries[j])
            end
        end
    end

    return pop_new, new_individuals
end

"""
    gadecimal_criteria(df::AbstractVector; ucb::Int, mean_tol::AbstractVector, std_tol::AbstractVector)

## Description
Criteria for the decimal-coded GA algorithm.
The criteria is based on the mean stability of the upper class boundary (3 last generations are considered) and on its standard deviation.
## Arguments
- df: vector with all the dataframes obtained during the process (check experiments.jl)
- ucb: upper class boundary
- mean_tol: vector with the tolerance for the mean values
- std_tol: vector with the tolerance for the std values
"""
function gadecimal_criteria(df::AbstractVector; ucb::Int, mean_tol::AbstractVector, std_tol::AbstractVector)

    n = length(df)
    i = n - 1
    k = n - 2

    if k <= 0
        return false
    else
        ncoord = ncol(df[end]) - 1
        mean_end = [mean(df[end][1:ucb,q]) for q in 1:ncoord]
        std_end = [std(df[end][1:ucb,q]) for q in 1:ncoord]
        mean_i = [mean(df[i][1:ucb,q]) for q in 1:ncoord]
        mean_k = [mean(df[k][1:ucb,q]) for q in 1:ncoord]

        Δmean_i = abs.(mean_end - mean_i)
        Δmean_k = abs.(mean_end - mean_k)

        if all(Δmean_i .< mean_tol) && all(Δmean_k .< mean_tol) && all(std_end .< std_tol)
            return true
        else
            return false
        end
    end
end
