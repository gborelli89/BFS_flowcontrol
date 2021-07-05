# -----------------------------------------------------------------------
# Codes used during the experiments
# -----------------------------------------------------------------------
# Required packages
using PyCall, Dates, JLD
include("cross_entropy.jl")
include("ga_decimal.jl")

# Including path to scanivalve 
push!(pyimport("sys")["path"], "path to scanivalve folder")


#---------------------------------------------------------------------------------------------------------------
# SCANIVALVE CONNECTION AND CONFIGURING
#---------------------------------------------------------------------------------------------------------------

"""
    scani_open(;ip="191.30.80.131")

## Description
Open connection with Scanivalve DSA 3217
## Arguments
- ip: IP device
## Example
```jldoctest
julia> dev = open_scani();

# Example of configuration
julia> dev.config(AVG=1, PERIOD=1000, FPS=1000) # for 16.0 s acquisition
```
"""
function scani_open(;ip="191.30.80.131")
    scani = pyimport("scanivalve")
    s = scani.Scanivalve(ip=ip)
    return s
end

"""
    scani_acquire(s; sname=missing)

## Description
Function to acquire and save measured data
## Arguments
- s: scanivalve connection already configured
- sname: name of the file to be saved. If missing, then only returns the pressure, the frequency and the DateTime of the measurements
## Example
```jldoctest
julia> s = scani_open()

julia> s.config(AVG=1, PERIOD=1000, FPS=1000)

julia> scani_acquire(s, "teste.jld")
```
"""
function scani_acquire(s; sname=missing)
    p, f = s.acquire()
    h = now()
    
    if ismissing(sname)
        return p, f, h
    else
        save(sname, "pressure", p, "freq", f, "h", h)
        return p, f, h
    end
end    


#---------------------------------------------------------------------------------------------------------------
# AUTOMATIC TEST WITH PREDEFINED positions
#---------------------------------------------------------------------------------------------------------------

"""
    testmatrix_onecylinder(;xlims::AbstractVector, ylim::AbstractVector, dx, dy)

## Description
Create a matrix with the measurements points for a test with one cylinder (no GA)
## Arguments
- xlims: vector with the low and the up boundaries for x direction
- ylims: vector with the low and the up boundaries for y direction
- dx, dy: x and y steps
## Example
```jldoctest
julia> testmatrix_onecylinder(xlims=[0,5],ylims=[-1,1], dx=5, dy=2)
4-element Array{OneCylType,1}:
 OneCylType(0.0, -1.0)
 OneCylType(0.0, 1.0)
 OneCylType(5.0, -1.0)
 OneCylType(5.0, 1.0)
```
"""
function testmatrix_onecylinder(;xlims::AbstractVector, ylims::AbstractVector, dx, dy)

    x = range(xlims[1], stop=xlims[2], step=dx)
    y = range(ylims[1], stop=ylims[2], step=dy)

    p = [OneCylType(i,j) for i in x for j in y]

    return p
end

"""
    autotest_ocylinder(s, robot, pos_vec::Vector{OneCylType}; sname="testdata.jld")

## Description
Perform multiple measurements in a given grid of positions (for 1 cylinder only).
## Arguments
- s: scanivalve connection configured
- robot: robot connection configured
- pos_vec: vector of OneCylType - grid of positions
- returntozero: if true then the cylinder returns back to the origin before new position
- sname: file name to be saved
## Example
```jldoctest
julia> s = scani_open();

julia> s.config(AVG=1, PERIOD=1000, FPS=1000);

julia> dev = stepper_init(1)

julia> grid = testmatrix_onecylinder(xlims=[0.0,90.0], ylims=[-30.0,20.0], dx=15.0, dy=5.0);

julia> d,f = autotest_1cylinder(s, dev, grid)
```
"""
function autotest_onecylinder(s, robot, pos_vec::Vector{OneCylType}; returntozero=false, sname="testdata.jld")

    coords = convert.(Array,pos_vec)
    pressure = Array{Array{Float64,2},1}()
    freq = Array{Float64,1}()
    testdate = Array{DateTime,1}()

    for pos in pos_vec
        println("Moving to position: x = ", pos.x, ", y = ", pos.y, "\n")
        if returntozero
            movecyl!(robot, OneCylType(0.0,0.0))
        end
        movecyl!(robot, pos)
        println("Acquiring data...\n")
        p, f, h = scani_acquire(s)
        push!(pressure, p)
        push!(freq, f)
        push!(testdate, h)
        println("Saving on file...\n")
        save(sname, "coords", coords, "pressure", pressure, "freq", freq, "h", testdate)    
    end

    println("End of test")

    return pressure, freq, testdate
end

#---------------------------------------------------------------------------------------------------------------
# CROSS ENTROPY METHOD 
#---------------------------------------------------------------------------------------------------------------

"""
    crossentropy(s, robot, pos_vec::AbstractVector, nbest::Int; mean_tol::AbstractVector, std_tol::AbstractVector, maxiter=10, dx=5.0, dy=0.0)

## Description
Cross entropy test for one cylinder
## Arguments
- sname: name of the file to be saved
- s: scanivalve connection configured
- robot: robot connection configured
- N: number of positions for the second iteration and on (for the first operation pos_vec is used)
- pos_vec: initial positions vector
- nbest: number of values which should be selected for the CE
- mean_tol: tolerance criteria for mean variation (for each direction)
- std_tol: tolerance criteria for standard deviation (for each direction)
- maxiter: maximum number of iterations 
- dx, dy: horizontal and vertical gap between cylinders at the origin
## Output
Array with all the positions and mean step base pressure as well as the best position found
"""
function crossentropy(sname, s, robot, N::Int, pos_vec::AbstractVector, nbest::Int; mean_tol::AbstractVector, 
                        std_tol::AbstractVector, maxiter=10)

    iter = 1
    df = []
    data = convert(DataFrame, pos_vec)
    data[!,:pm] = repeat([0.0], nrow(data))

    best = DataFrame()

    while iter <= maxiter
        println("Iteration ", iter)
        pm = []

        for pos in pos_vec
            println("Moving to position ", convert(Array, pos), "\n")
            strategy = beststrategy(robot, pos)
            movecyl!(robot, pos, strategy=strategy)
            println("Acquiring data...\n")
            p,_,_ = scani_acquire(s)
            push!(pm, mean(p[:,1:9]))
        end
        
        data.pm = pm
        push!(df, deepcopy(data))
        append!(best, data)
        best = best[best.pm .== maximum(best.pm),:]
                                                                                
        println(data)

        lims, converge = crossentropy(df, nbest; mean_tol=mean_tol, std_tol=std_tol)
        println(lims)
        println("converged? ", converge)

        if converge
            return df, best
            break
        else
            pos_vec = randpos(N, lims)
            data = convert(DataFrame, pos_vec)
            data[!,:pm] = repeat([0.0], nrow(data))        
        end

        save(sname, "df", df, "best", best)
        iter += 1
    end

    return df, best
end


#---------------------------------------------------------------------------------------------------------------
# GENETIC ALGORITHM 
#---------------------------------------------------------------------------------------------------------------
#
# DECIMAL-CODED
#
"""
    gadecimal(s, robot, pop_init::AbstractVector; Tm::AbstractVector, α0::AbstractVector, maxgen::Int, ucb::Int, lcb::Int, digits=0, boundaries::AbstractVector, mean_tol=AbstractVector, std_tol=AbstractVector, dx=5.0, dy=0.0, measureall=false)

## Description 
Decimal-coded, mutation based genetic ALGORITHM
## Arguments
- s: scanivalve connection configured
- robot: robot connection configured
- pos_vec: initial position vector (CylType)
- Tm: mutation rate for each direction
- α0: shift for the first generation
- maxgen: total number of generations
- ucb: upper class boundary (elitism)
- lcb: lower class boundary
- digits: decimal places
- boundaries: domain limits for the direction
- mean_tol: vector with the tolerance for the mean values
- std_tol: vector with the tolerance for the std values
- dx, dy: horizontal and vertical gap between cylinders at the origin
- measureall: if true then all the positions are measured each iteration (even repeated ones)
"""
function gadecimal(sname, s, robot, pos_vec::AbstractVector; Tm::AbstractVector, α0::AbstractVector, maxgen::Int, ucb::Int, 
        lcb::Int, digits=0, boundaries::AbstractVector, mean_tol::AbstractVector, std_tol::AbstractVector, measureall=false)

    gen = 1
    N = length(pos_vec)
    df = []
    tomeasure = 1:N

    df = []
    data = convert(DataFrame, pos_vec)
    data[!,:fitness] = repeat([0.0], N)

    while gen <= maxgen
        println("Generation ", gen)
        pm = data.fitness
        
        for k in tomeasure
            println("Moving to position ", convert(Array, pos_vec[k]), "\n")
            strategy = beststrategy(robot, pos_vec[k])
            movecyl!(robot, pos_vec[k], strategy=strategy)
            println("Acquiring data...\n")
            p,_,_ = scani_acquire(s)
            pm[k] = mean(p[:,1:9])
        end
        
        data.fitness = pm
        push!(df, deepcopy(data))

        println(data)

        converge = gadecimal_criteria(df, ucb=ucb, mean_tol=mean_tol, std_tol=std_tol)
        if converge
            return df
            break
        end

        data,tomeasure = gadecimal_step(data, Tm=Tm, α0=α0, gen=gen, maxgen=maxgen, ucb=ucb, lcb=lcb, 
                        digits=digits, boundaries=boundaries)

        if measureall
            tomeasure = 1:N
        end

        pos_vec = convert(eltype(pos_vec), data[!,Not(:fitness)])
        
        save(sname, "df", df)
        gen += 1
    end

    return df
end
