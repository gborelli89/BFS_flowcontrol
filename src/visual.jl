# -----------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------
include("cylinder_structure.jl")
using Plots
gr()

"""
    circleShape(xc, yc, r)

## Description
Create a circunference
## Arguments
- xc and yc: center
- r: radius
"""
function circleShape(xc, yc, r)
    θ = LinRange(0,2*π, 500)
    @. xc + r*sin(θ), yc + r*cos(θ)
end 

"""
    function plotCyl(pos::Array{OneCylType,1}; cylAlpha=0.5, B=(-50,20), L=(-10,100), r=1.2, dmfac=50.0, δc=(4,0), ptitle="")

## Description
Plot one cylinder (multiple position)
## Arguments
- pos: array of OneCylType
- cylAlpha: transparency paramenter
- B: tuple with x limits
- L: tuple with y limits
- r: cylinder radius
- dmfac: dimensionless factor
- δc: origin distance from the step
- ptitle: plot title
"""
function plotCyl(pos::Array{OneCylType,1}; cylAlpha=0.5, B=(-50,20), L=(-10,100), r=1.2, 
    dmfac=50.0, δc=(4,0), ptitle="")
   
    n = length(pos)
    xb = [L[1],0.0,0.0] ./ dmfac
    yb = [0.0,0.0,B[1]] ./ dmfac
    r = r / dmfac
    δc = δc ./ dmfac
    L = L ./ dmfac
    B = B ./ dmfac

    if length(cylAlpha) == 1
        cylAlpha = repeat([cylAlpha], n)
    end

    p = plot(xb,yb, lw=2, xlims=L, ylims=B, c=:gray, legend=false, aspect_ratio=1)
    xlabel!("x/B")
    ylabel!("y/B")
    title!(ptitle)
   
    for i in 1:n
        plot!(circleShape(pos[i].x/dmfac + δc[1], pos[i].y/dmfac + δc[2], r), seriestype=[:shape], lw=0.1, 
            c=:blue, linecolor=false, fillalpha=cylAlpha[i])
    end
   
    display(p)
end

"""
    function plotCyl(pos::Array{TwoCylType,1}; cylAlpha=[0.5], B=(-50,20), L=(-10,120), r=1.2, dmfac=50.0, δc=(4,0,10,0), ptitle="")

## Description
Plot two cylinder (multiple position)
## Arguments
- pos: array of TwoCylType
- cylAlpha: transparency paramenter
- B: tuple with x limits
- L: tuple with y limits
- r: cylinder radius
- dmfac: dimensionless factor
- δc: array with the origin distance from the step and the distance between cylinders at the origin
- ptitle: plot title
"""
function plotCyl(pos::Array{TwoCylType,1}; cylAlpha=[0.5], B=(-50,20), L=(-10,120), r=1.2, 
    dmfac=50.0, δc=(4,0,10,0), ptitle="")
   
    n = length(pos)
    xb = [L[1],0.0,0.0] ./ dmfac
    yb = [0.0,0.0,B[1]] ./ dmfac
    r = r / dmfac
    δc = δc ./ dmfac
    B = B ./ dmfac
    L = L ./ dmfac

    p = plot(xb,yb, lw=2, xlims=L, ylims=B, c=:gray, legend=false, aspect_ratio=1)
    xlabel!("x/B")
    ylabel!("y/B")
    title!(ptitle)

    if length(cylAlpha) == 1
        cylAlpha = repeat([cylAlpha], n)
    end

    for i in 1:n
        Px = [pos[i].x/dmfac + δc[1], pos[i].x/dmfac + pos[i].dx/dmfac + δc[3]]
        Py = [pos[i].y/dmfac + δc[2], pos[i].y/dmfac + pos[i].dy/dmfac + δc[4]]
        plot!(circleShape(Px[1], Py[1], r), seriestype=[:shape], lw=0.1, 
            c=:blue, linecolor=false, fillalpha=cylAlpha[i])
        plot!(circleShape(Px[2], Py[2], r), seriestype=[:shape], lw=0.1, 
            c=:blue, linecolor=false, fillalpha=cylAlpha[i])
        plot!(Px, Py, c=:red)
    end

    display(p)
end 

"""
    function plotCyl(pos::Array{ThreeCylType,1}; cylAlpha=[0.5], B=(-50,20), L=(-10,130), r=1.2, dmfac= 50.0, δc=(4,0,10,0,16,0), ptitle="")

## Description
Plot three cylinder (multiple position)
## Arguments
- pos: array of ThreeCylType
- cylAlpha: transparency paramenter
- B: tuple with x limits
- L: tuple with y limits
- r: cylinder radius
- dmfac: dimensionless factor
- δc: array with the origin distance from the step and the distance between cylinders at the origin
- ptitle: plot title
"""
function plotCyl(pos::Array{ThreeCylType,1}; cylAlpha=[0.5], B=(-50,20), L=(-10,130), r=1.2, 
    dmfac= 50.0, δc=(4,0,10,0,16,0), ptitle="")
    
    n = length(pos)
    xb = [L[1],0.0,0.0] ./ dmfac
    yb = [0.0,0.0,B[1]] ./ dmfac
    r = r / dmfac
    δc = δc ./ dmfac
    B = B ./ dmfac
    L = L ./ dmfac    

    p = plot(xb,yb, lw=2, xlims=L, ylims=B, c=:gray, legend=false, aspect_ratio=1)
    xlabel!("x/B")
    ylabel!("y/B")
    title!(ptitle)
    
    if length(cylAlpha) == 1
        cylAlpha = repeat([cylAlpha], n)
    end

    for i in 1:n
        Px = [pos[i].x/dmfac + δc[1], pos[i].x/dmfac + pos[i].dx1/dmfac + δc[3], pos[i].x/dmfac + pos[i].dx1/dmfac + pos[i].dx2/dmfac + δc[5]]
        Py = [pos[i].y/dmfac + δc[2], pos[i].y/dmfac + pos[i].dy1/dmfac + δc[4], pos[i].y/dmfac + pos[i].dy1/dmfac + pos[i].dy2/dmfac + δc[6]]
        plot!(circleShape(Px[1], Py[1], r), seriestype=[:shape], lw=0.1, 
            c=:blue, linecolor=false, fillalpha=cylAlpha[i])
        plot!(circleShape(Px[2], Py[2], r), seriestype=[:shape], lw=0.1, 
            c=:blue, linecolor=false, fillalpha=cylAlpha[i])
        plot!(circleShape(Px[3], Py[3], r), seriestype=[:shape], lw=0.1, 
            c=:blue, linecolor=false, fillalpha=cylAlpha[i])
        plot!(Px, Py, c=:red)
    end

    display(p)
end 
# Methods for one position
function plotCyl(pos::CylType; cylAlpha=0.5, B=(-50,20), L=(-10,130), r=1.2)
    plotCyl([pos], cylAlpha=cylAlpha, B=B, L=L, r=r)
end

"""
    plotGeneration(df::DataFrame, rmcol=missing)

## Description
Plot cylinders for dataframe data format
## Arguments
- df: dataframe with x, y, dx, dy parameters
- rmcol: columns to be removed
"""
function plotGeneration(df::DataFrame; rmcol=missing, dmfac=50.0, ptitle="")

    if !ismissing(rmcol)
        df = df[!,Not(rmcol)]
    end

    nc = ncol(df)
    if nc == 2
        el = convert(OneCylType, df)
    elseif nc == 4
        el = convert(TwoCylType, df)
    elseif nc == 6
        el = convert(ThreeCylType, df)
    else
        throw(DimensionMismatch)
    end

    plotCyl(el, dmfac=dmfac, ptitle=ptitle)
end