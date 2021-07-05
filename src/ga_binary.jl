# -----------------------------------------------------------------------
# Binary coded genetic algorithm
# -----------------------------------------------------------------------
# Required packages
using BitBasis
using Statistics
include("stepper_move.jl")
include("optim_move.jl")
include("random_init.jl")
include("visual.jl")
include("simdata.jl") 
# -----------------------------------------------------------------------
# Important onbservations
# These codes were not actually used in the final experiments
# Instead, the decimal coded algorithm was preferred
# It might not work well. 
# The simulator implemented works for one cylinder onlys
# -----------------------------------------------------------------------
# Useful functions
# -----------------------------------------------------------------------

"""
    cylpos2bits(cyl::CylType; size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4])

## Description
Convert a position of any CylType to bit.
## Arguments
- cyl: position of type OneCylType, TwoCylType or ThreeCylType
- size: size of the search 2^size (mm)
- coord_zero: initial coordinates (offset)
## Example
```jldoctest
julia> cylpos2bits(ThreeCylType(10,2,5,0,6,-1))
6-element Array{BitArray{1},1}:
 [0, 1, 0, 1, 0]
 [1, 1, 1, 0]
 [1, 0, 1]
 [0, 0, 1]
 [0, 1, 1]
 [1, 1, 0]
```
"""
function cylpos2bits(cyl::CylType; size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4])
    a = convert(Array, cyl)
    n = length(a)
    id = 1:n

    # Applying offset on the coordinates 
    @. a[id] = a[id] - coord_zero[id]

    # convert to bit
    a_bit = bitarray.(a, size[id])
    
    return a_bit
end


"""
    bits2cylpos(b::AbstractArray; coord_zero=[0,-5,0,-4,0,-4])

## Description
Convert a position of any CylType to bit.
## Arguments
- b: position on bits
- coord_zero: initial coordinates (offset) 
## Example
```jldoctest
julia> b = cylpos2bits(ThreeCylType(10,2,5,0,6,-1));

julia> bits2cylpos(b)
ThreeCylType(10.0, 2.0, 5.0, 0.0, 6.0, -1.0)
```
"""
function bits2cylpos(b::AbstractArray; coord_zero=[0,-5,0,-4,0,-4])
    n = length(b)
    id = 1:n

    p = packbits.(b)

    # Applying offset on the coordinates 
    @. p[id] = p[id] + coord_zero[id]

    if n == 2
        pos = OneCylType(p[1],p[2])
    elseif n == 4
        pos = TwoCylType(p[1],p[2],p[3],p[4])
    elseif n == 6
        pos = ThreeCylType(p[1],p[2],p[3],p[4],p[5],p[6])
    else
        error("Number of cylinders not defined")
    end
    
    return pos
end

"""
    sigmoid(gen; offset, α)

## Description
Sigmoidal function to control the mutation.
## Arguments
- gen: number of generation
- offset: generation offset parameter (20)
- α: power factor (0.4)
"""
sigmoid(gen; offset, α) = 1 - 1/(1+exp(-α*(gen-offset)))


"""
    dict_init(type=ThreeCylType)

## Description
Initiate a dictionary for the read data.
## Arguments
- type: key type, use OneCylType, TwoCylType or ThreeCylType
"""
dict_init(type=ThreeCylType) = Dict{type, Array{Float64,1}}()

"""
    repeatedpos(dict, pos::CylType)

## Description
Find if the position was not measured. If measured, add one count on the dictionary.
## Arguments
- dict: dictionary with the positions and the data measured
- pos: position of CylType 
"""
function repeatedpos(dict::Dict, pos::CylType)
    if haskey(dict, pos)
        dict[pos][2] += 1
        return false
    else
        return true
    end
end

"""
    newposition(dict::Dict, posArray::ArrayCylType)

## Description
Discover which positions of the generation must be measured.
## Arguments
- dict: dictionary with the positions and the data measured
- posArray: array with the positions of CylType
"""
function newposition(dict::Dict, posArray::ArrayCylType)
    id = [repeatedpos(dict,i) for i in posArray]
    new_posArray = posArray[id]
    return new_posArray
end

# -----------------------------------------------------------------------
# Convergence
# -----------------------------------------------------------------------

"""
    criteriarule(new_coord::AbstractVector, old_coord::AbstractVector; μ, σ)

## Description
General criteria rule.
## Arguments
- old_coord: array with floats of the old coordinate
- new_coord: array with floats of the new coordinate
- μ: mean parameter
- σ: standard deviation parameter
"""
function criteriarule(new_coord::AbstractVector, old_coord::AbstractVector; μ, σ)
    
    μ = float(μ)
    σ = float(σ)

    Δcoord = abs(mean(new_coord)-mean(old_coord))

    if Δcoord < μ && std(new_coord) < σ 
        return true
    else
        return false
    end
end

"""
    criteria(new_gen::Array{T,1}, old_gen::Array{T,1}; μ, σ, Δμ, Δσ) where T<:CylType

## Description
Criteria for one, two or three control cylinders.
## Arguments
- old_gen: last generation measured
- new_gen: new generation measured
- μ, Δμ: mean parameters
- σ, Δσ: standard deviation parameters
"""
function criteria(new_gen::Array{OneCylType,1}, old_gen::Array{OneCylType,1}; μ, σ, Δμ=0.0, Δσ=0.0)
   
    n = length(new_gen)
   
    x_old = [old_gen[i].x for i in 1:n]
    x_new = [new_gen[i].x for i in 1:n]
   
    if criteriarule(x_new, x_old , μ=μ, σ=σ)
        y_old = [old_gen[i].y for i in 1:n]
        y_new = [new_gen[i].y for i in 1:n]
        return criteriarule(y_new, y_old, μ=μ, σ=σ)
    else
        return false
    end
end

function criteria(new_gen::Array{TwoCylType,1}, old_gen::Array{TwoCylType,1}; μ, σ, Δμ, Δσ)
   
    n = length(new_gen)
   
    x_old = [old_gen[i].x for i in 1:n]
    x_new = [new_gen[i].x for i in 1:n]
   
    if criteriarule(x_new, x_old , μ=μ, σ=σ)
        y_old = [old_gen[i].y for i in 1:n]
        y_new = [new_gen[i].y for i in 1:n]
        
        if criteriarule(y_new, y_old, μ=μ, σ=σ)
            dx_old = [old_gen[i].dx for i in 1:n]
            dx_new = [new_gen[i].dx for i in 1:n]

            if criteriarule(dx_new, dx_old, μ=Δμ, σ=Δσ)
                dy_old = [old_gen[i].dy for i in 1:n]
                dy_new = [new_gen[i].dy for i in 1:n]
                return criteriarule(dy_new, dy_old, μ=Δμ, σ=Δσ)
            else 
                return false
            end
        else
            return false
        end
    else
        return false
    end
end

function criteria(new_gen::Array{ThreeCylType,1}, old_gen::Array{ThreeCylType,1}; μ, σ, Δμ, Δσ)
   
    n = length(new_gen)
   
    x_old = [old_gen[i].x for i in 1:n]
    x_new = [new_gen[i].x for i in 1:n]
   
    if criteriarule(x_new, x_old , μ=μ, σ=σ)
        y_old = [old_gen[i].y for i in 1:n]
        y_new = [new_gen[i].y for i in 1:n]
        
        if criteriarule(y_new, y_old, μ=μ, σ=σ)
            dx1_old = [old_gen[i].dx1 for i in 1:n]
            dx1_new = [new_gen[i].dx1 for i in 1:n]

            if criteriarule(dx1_new, dx1_old, μ=Δμ, σ=Δσ)
                dy1_old = [old_gen[i].dy1 for i in 1:n]
                dy1_new = [new_gen[i].dy1 for i in 1:n]
                
                if criteriarule(dy1_new, dy1_old, μ=Δμ, σ=Δσ)
                    dx2_old = [old_gen[i].dx2 for i in 1:n]
                    dx2_new = [new_gen[i].dx2 for i in 1:n]

                    if criteriarule(dx2_new, dx2_old, μ=Δμ, σ=Δσ)
                        dy2_old = [old_gen[i].dy2 for i in 1:n]
                        dy2_new = [new_gen[i].dy2 for i in 1:n]
                        return criteriarule(dy2_new, dy2_old, μ=Δμ, σ=Δσ)
                    else
                        return false
                    end
                else
                    return false
                end
            else 
                return false
            end
        else
            return false
        end
    else
        return false
    end
end

# -----------------------------------------------------------------------
# Genetic Algorithm functions
# -----------------------------------------------------------------------

"""
    remove_individual(rm, pop)

## Description
Remove individual from the population.
## Arguments
- rm: individual removed
- pop: population
"""
function remove_individual(rm, pop)
    id = [i != rm for i in pop]
    return pop[id]
end

"""
    roulette_selection(val, pop)

## Description
Selection of the roulette type.
## Arguments
- val: evaluation of each individual, according to a feature
- pop: population
"""
function roulette_selection(val, pop)
    fitness = round.(cumsum(val/sum(val)), digits=5)
    n = length(val) - 1 # one individual already chosen (eletive selection)
    id = []

    for i in 1:n
        idtemp = sum(rand() .> fitness) + 1
        id = push!(id, idtemp)
    end

    return pop[id]
end

"""
    crossover(parent1::Array{BitArray{1},1}, parent2::Array{BitArray{1},1}; size)

## Description
Crossover.
## Arguments
- parent1: first individual to crossover (bit)
- parent2: second individual to crossover (bit)
- size: bit sizes
"""
function crossover(parent1::Array{BitArray{1},1}, parent2::Array{BitArray{1},1}; size)
    
    ncoord = length(parent1)
    cut = vcat(rand([1,2,3], 2), rand([1,2],ncoord-2)) # can be changed
    child1 = Array{BitArray{1},1}()
    child2 = Array{BitArray{1},1}()

    for i in 1:ncoord
        k = cut[i] + 1
        push!(child1, vcat(parent1[i][1:cut[i]], parent2[i][k:size[i]]))
        push!(child2, vcat(parent2[i][1:cut[i]], parent1[i][k:size[i]]))
    end

    return [child1, child2]

end

"""
    mutation(ind::Array{BitArray{1},1}; Tm)

## Description
Mutation.
## Arguments
- ind: single individual of the population
- Tm: mutation rate
"""
function mutation(ind::Array{BitArray{1},1}; Tm)
    
    for i in 1:length(ind)
        idm = rand(length(ind[i])) .< Tm
        ind[i][idm] = .!ind[i][idm]
    end

    return ind
end


"""
    ga_onestep(pop, dict; Tr=0.25, Tm=0.005, size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4])

## Description
Implementation of the genetic algorithm (one iteration).
## Arguments
- pop: population
- Tr: crossover rate
- Tm: mutation rate
- size: size of the search 2^size (mm)
- coord_zero: initial coordinates (offset)
"""
function ga_onestep(pop, dict; Tr=0.25, Tm=0.005, size=[5,4,3,3,3,3], coord_zero=[0,-5,0,-4,0,-4])

    # Fitness of each individual
    val = [dict[i][1] for i in pop]
    #println("fitness ok")

    # Eletive selection
    idmax = argmax(val)
    elite = pop[idmax]
    #println("eletive ok")

    # Roulette selection
    pop_sel = roulette_selection(val, pop)
    #println("roulette ok")

    # Convert to bit
    pop_bit = cylpos2bits.(pop_sel, size=size, coord_zero=coord_zero)
    np_bit = length(pop_bit)
    #println("convert to bit ok")

    # Crossover
    idr = rand(np_bit) .< Tr
    idr_inv = .!(idr)
    
    if sum(idr) >= 2
        parents = pop_bit[idr]
        pop_bit_r = pop_bit[idr_inv]
        
        for i in 1:2:sum(idr)
            k = i+1
            if i == sum(idr)
                push!(pop_bit_r, parents[i])
                break
            end
            child1, child2 = crossover(parents[i], parents[k], size=size)
            push!(pop_bit_r, child1, child2)
        end
        
    else
        pop_bit_r = pop_bit
    end
    #println("crossover ok")

    # Mutation 
    pop_bit_m = mutation.(pop_bit_r, Tm=Tm)
    #println("mutation ok")

    # New generation
    new_gen = bits2cylpos.(pop_bit_m, coord_zero=coord_zero)
    push!(new_gen, elite) # adding eletive selection

    #return [elite, pop_sel, pop_bit, pop_bit_r, pop_bit_m, new_gen]
    return new_gen

end


"""
    ga_simulator(;n=30, nc=3, Tr=0.25, Tmmax=0.01, offset=25, α=0.4, μ=1.0, σ=0.5, Δμ=1.0, Δσ=1.0, ngen_max=100)

## Description
GA simulator.
## Arguments
- n: number of individuals in population
- nc: number of cylinders
- Tr: crossover rate
- Tmmax: maximum mutation rate
- offset and α: sigmoid parameters (for mutation decay)
- μ, σ, Δμ, Δσ: convergence criteria parameters
- ngen_max: max number of generations allowed
"""
function ga_simulator(;n=30, nc=3, Tr=0.25, Tmmax=0.01, offset=25, α=0.4, μ=1.0, σ=0.5, Δμ=1.0, Δσ=1.0, ngen_max=100)
    
    size = [5,4,3,3,3,3]
    coord_zero=[0,-5,0,-4,0,-4]
    pop = population_init(n, nc=nc, size=size[1:2*nc], coord_zero=coord_zero[1:2*nc]) # initializing population
    dict = dict_init(typeof(pop[1])) # initializing dictionary

    data = simDataOne("data_test01") # get simulation data

    ngen = 1 # initializing generation number
    criteriaConv = false

    new_pop = pop
    NI = 0 # total number of individuals (visited at least once)
    NS = 0 # total number of steps
    TT = 0 # total time
    pos = randindividual(nc, size=size) # initial position of the robot

    while ngen < ngen_max && !criteriaConv

        println("Geração: ", ngen)

        toMeasure = newposition(dict, new_pop) # find which individuals must be evaluated
        NI += length(toMeasure)
        
        if length(toMeasure) > 0
            if length(toMeasure) == 1
                toMeasure_opt = toMeasure
            else
                toMeasure_opt = toMeasure[greedypath(pos, toMeasure)]
            end 
            
            #toMeasure_opt = toMeasure # not applying path optimization

            for i in 1:length(toMeasure_opt)
                NS += sum(abs.(coord2steps(pos,toMeasure_opt[i])[1].step)) # find total steps for each individual
                pos = toMeasure_opt[i]
            end

            plotCyl_Gen(toMeasure_opt; B=(-20,20), L=(-20,70), r=1.2) # plot position measured
            [dict[i] = [data[i], 1] for i in toMeasure_opt] # collecting simulation data
        end

        TT = 2*1.29*nc*NI + 0.002*NS + 15*NI # find total time (seconds)
        println("    new individuals: ", NI)
        println("    total steps (stepper motor): ", NS)
        println("    simulation time: ", round(TT/60,digits=2), " min")

        old_pop = new_pop
        Tm = Tmmax*sigmoid(ngen, offset=offset, α=α) # finding mutation rate
        new_pop = ga_onestep(old_pop, dict, Tr=Tr, Tm=Tm, size=size[1:2*nc], coord_zero=coord_zero[1:2*nc]) # run genetic algorithm

        criteriaConv = criteria(new_pop, old_pop, μ=μ, σ=σ, Δμ=Δμ, Δσ=Δσ) # evaluate criteria
        ngen += 1 # move to next generation
        println("    criteria: ", criteriaConv)

    end

    return [new_pop, dict, TT, NI, NS, criteriaConv]
end


# Determining all the parameters for one cylinder case study
# ga_simulator(;n=30, nc=3, Tr=0.25, Tmmax=0.01, offset=25, α=0.4, 
#    μ=1.0, σ=0.5, Δμ=1.0, Δσ=1.0, ngen_max=100)

function params_onecylinder(ncases=100; maxval= -62.74781119, nrep = 200, ngen_max=1000)
    
    Tmean = []; Tsd = []; V = []; C = []; params = []
    
    for i in 1:ncases
        pop_size = rand(10:30)
        tr_range = rand(0.1:0.05:0.6)
        tm_range = rand(0.002:0.002:0.03)
        sigmoid_offset = rand(10:50)
        sigmoid_alpha = rand(0.1:0.05:0.9)
        convergence_mu = rand(0.2:0.1:2.0)
        convergence_sigma = rand(0.1:0.05:0.9)

        push!(params, [pop_size,tr_range,tm_range,sigmoid_offset,sigmoid_alpha,convergence_mu,convergence_sigma])
        
        println("Iteration ", i)
        
        Ttemp = []; Vtemp = []
        nconv = 0
        for _ in 1:nrep
            _,d,t,conv = ga_simulator(n=pop_size, nc=1, Tr=tr_range, Tmmax=tm_range, offset=sigmoid_alpha,
                        α=sigmoid_alpha, μ=convergence_mu, σ=convergence_sigma, ngen_max=ngen_max)
            
            if conv 
                push!(Ttemp, t/60)
                v = d[argmax(d)][1]
                diff_v = abs(v-maxval)
                push!(Vtemp, diff_v)
                nconv += 1
            end
    
        end

        println("   End of iteration")
        println("   Summary -> total: ", nrep, ", converged: ", nconv)

        push!(Tmean, mean(Ttemp))
        push!(Tsd, std(Ttemp))
        push!(V, mean(Vtemp))
        push!(C, nconv)
    end
    
    return [Tmean,Tsd,V,C,params]
end






    
    
