# -----------------------------------------------------------------------
# Functions that provide the control of alls the stepper motors
# ----------------------------------------------------------------------- 
# Required packages
using StepperControl, Optim
using BSON: @load
include("cylinder_structure.jl")

"""
    step2coord_x(;sprx=2048, R=80.0, xorig=-41.5)

## Description
Generates a  function which converts steps into displacement in x direction
## Arguments
- sprx: steps per revolution in x system
- R: bell crank radius
- xorig: horizontal distance from the center of rotation to the origin
## Example
```jldoctest
julia> f = step2coord_x();

julia> f(20)
4.323700717312371
```
"""
function step2coord_x(;sprx=2048, R=80.0, xorig=-41.5)
    θ0 = asin(xorig/R)
    x(v) = R*(sin(θ0 + 2*π*v/sprx) - sin(θ0))
    return x
end


"""
    coord2step_x(;sprx=2048, R=80.0, xorig=-41.5)

## Description
Generates a  function which converts displacement in x direction into steps 
## Arguments
- sprx: steps per revolution in x system
- R: bell crank radius
- xorig: horizontal distance from the center of rotation to the origin
## Example
```jldoctest
julia> f = coord2step_x();

julia> f(4.32)
20
```
"""
function coord2step_x(;sprx=2048, R=80.0, xorig=-41.5)
    θ0 = asin(xorig/R)
    x(v) = Int(round( (sprx/(2π)) * (asin(v/R + sin(θ0)) - θ0) ))
    return x
end


"""
    step2coord_y(;sprx=2048, spry=2048, R=80.0, r=12.5, xorig=-41.5)

## Description
Generates a  function which converts steps into displacement in y direction
## Arguments
- sprx: steps per revolution in x system
- spry: steps per revolution in y system
- R: bell crank radius
- r: gear radius
- xorig: horizontal distance from the center of rotation to the origin
## Example
```jldoctest
julia> f = step2coord_y();

julia> f([20,100])
1.5861345582862714
```
"""
function step2coord_y(;sprx=2048, spry=2048, R=80.0, r=12.5, xorig=-41.5)
    θ0 = asin(xorig/R)
    y(v) = 2π*v[2]*r/spry - R*(cos(θ0 + 2*π*v[1]/sprx) - cos(θ0))
    return y
end


"""
    coord2step_y(;sprx=2048, spry=2048, R=80.0, r=12.5, xorig=-41.5)

## Description
Generates a  function which converts displacement in y direction into steps 
## Arguments
- sprx: steps per revolution in x system
- spry: steps per revolution in y system
- R: bell crank radius
- r: gear radius
- xorig: horizontal distance from the center of rotation to the origin
## Example
```jldoctest
julia> f = coord2step_y();

julia> f(4.32, 1.58)
100
```
"""
function coord2step_y(;sprx=2048, spry=2048, R=80.0, r=12.5, xorig=-41.5)
    θ0 = asin(xorig/R)
    y(v) = Int(round( (spry/(2π*r)) * (v[2] + R*(cos(asin(v[1]/R + sin(θ0))) - cos(θ0))) ))
    return y
end


"""
    beststrategy(dev::StepperSystem, posnew::OneCylType)

## Description
Find best moving strategy to avoid collisions
## Arguments
- dev: stepper system from StepperControl library
- posnew: new CylType position
"""
beststrategy(dev::StepperSystem, posnew::OneCylType) = 1:2

function beststrategy(dev::StepperSystem, posnew::TwoCylType)

    posold = getpos(dev)

    if posnew[1] >= posold[1]
        if posnew[1] >= 55.0
            strategy = [3,4,2,1]
        else
            strategy = [3,4,1,2]
        end
    else
        if posnew[1] <= 41.0
            strategy = [1,2,4,3]
        else
            strategy = 1:4
        end
    end
   
    return strategy
end 

function beststrategy(dev::StepperSystem, posnew::ThreeCylType)

    posold = getpos(dev)

    if posnew[1] >= posold[1]
        if posnew[1]+posnew[5] <= 75.0
            strategy = [5,6,3,4,2,1]
        else
            strategy = [5,6,4,3,2,1]
        end
    else
        if posnew[1] >= 42.0
            strategy = [1,2,3,4,6,5]
        else
            strategy = [1,2,4,3,6,5]
        end
    end
   
    return strategy
end 


"""
    open_onecylinder(sprx=2121, spry=2018; R=80.0, r=12.5, xorig=-41.5, testnocon=false)

## Description 
Function to get connection and configure stepper system for one cylinder
## Arguments
- sprx: steps per revolution in x direction
- spry: steps per revolution in y direction
- R: bell crank length
- r: gear radius
- xorig: horizontal distance from de center of the rotor to the origin
- testnocon: if true then no connection is given (for testing purposes)
"""
function open_onecylinder(sprx=2099, spry=2007; R=80.0, r=12.5, xorig=-41.5, testnocon=false)

    # getting connection
    if testnocon
        dev = stepper_open(2, testnocon=true)
    else
        dev = stepper_open(2)
    end 

    # transformation functions
    s2cx = step2coord_x(sprx=sprx, R=R, xorig=xorig)
    s2cy = step2coord_y(sprx=sprx, spry=spry, R=R, r=r, xorig=xorig)
    c2sx = coord2step_x(sprx=sprx, R=R, xorig=xorig)
    c2sy = coord2step_y(sprx=sprx, spry=spry, R=R, r=r, xorig=xorig)

    stepper_config!(dev, motorID=["x1","y1"], step2coord=[s2cx,s2cy], coord2step=[c2sx,c2sy], depend=[1,[1,2]])

    return dev
end


"""
    open_cylinderNN(NNmodelfile::AbstractString; xorig=-59.0, motorID=["x2","y2"], testnocon=false)

## Description 
Function to get connection and configure stepper system for one cylinder considering neural net models
## Arguments
- NN: neural network
- xorig: horizontal distance from de center of the rotor to the origin
- motorID: stepper motor ids
- testnocon: if true then no connection is given (for testing purposes)
"""
function open_cylinderNN(NN; xorig=-59.0, motorID=["x2","y2"], testnocon=false)

    # getting connection
    if testnocon
        dev = stepper_open(2, testnocon=true)
    else
        dev = stepper_open(2)
    end 
    
    # ideal models
    θ0 = asin(xorig/80)
    s2cx(sx) = 80*(sin(θ0+ 2π*sx/2094) - sin(θ0))
    s2cy(s) = 2π*12.5*s[2]/1994 - 80*(cos(θ0 + 2*π*s[1]/2094) - cos(θ0)) + NN(s[1])
    c2sx(x) = Int(round( (2094/(2π))*(asin(x/80+sin(θ0)) - θ0) ))
    
    approx_c2sy(c) = (1994/(2π*12.5))*(c[2] + 80*(cos(asin(c[1]/80+sin(θ0)))-cos(θ0)))
    function c2sy(c)
        sx = c2sx(c[1])
        minval(sy) = sum(abs2, c[2] - s2cy([sx,sy]))
        r = optimize(v -> minval(first(v)), [approx_c2sy(c)], LBFGS())
        return Int(round(r.minimizer[1]))
    end
        
    stepper_config!(dev, motorID=motorID, step2coord=[s2cx,s2cy], coord2step=[c2sx,c2sy], depend=[1,[1,2]])

    return dev
end

"""
    open_twocylinder(; testnocon=false)

"""
function open_twocylinder(; testnocon=false)

    # getting connection
    if testnocon
        dev = stepper_open(4, testnocon=true)
    else
        dev = stepper_open(4)
    end 

    # transformation functions cylinder 1
    θ0₁ = asin(-41.5/80.0)
    s2cx1(sx) = 80*(sin(θ0₁+ 2π*sx/2099) - sin(θ0₁))
    s2cy1(s) = 2π*12.5*s[2]/2007 - 80*(cos(θ0₁ + 2*π*s[1]/2099) - cos(θ0₁))
    c2sx1(x) = Int(round( (2099/(2π))*(asin(x/80+sin(θ0₁)) - θ0₁) ))    
    c2sy1(c) = Int(round( (2007/(2π*12.5))*(c[2] + 80*(cos(asin(c[1]/80+sin(θ0₁)))-cos(θ0₁))) ))

    # transformation functions cylinder 2
    θ0₂ = asin(-59.0/80.0)
    s2cx2(sx) = 80*(sin(θ0₂+ 2π*sx/2094) - sin(θ0₂))
    s2cy2(s) = 2π*12.5*s[2]/1994 - 80*(cos(θ0₂ + 2*π*s[1]/2094) - cos(θ0₂)) + NNyx2(s[1])
    c2sx2(x) = Int(round( (2094/(2π))*(asin(x/80+sin(θ0₂)) - θ0₂) ))    
    approx_c2sy2(c) = (1994/(2π*12.5))*(c[2] + 80*(cos(asin(c[1]/80+sin(θ0₂)))-cos(θ0₂)))
    function c2sy2(c)
        sx = c2sx2(c[1])
        minval(sy) = sum(abs2, c[2] - s2cy2([sx,sy]))
        r = optimize(v -> minval(first(v)), [approx_c2sy2(c)], LBFGS())
        return Int(round(r.minimizer[1]))
    end


    stepper_config!(dev, motorID=["x1","y1","x2","y2"], step2coord=[s2cx1,s2cy1,s2cx2,s2cy2], 
                    coord2step=[c2sx1,c2sy1,c2sx2,c2sy2], depend=[1,[1,2],3,[3,4]])

    return dev
end


"""
    open_threecylinder(; testnocon=false)
    
"""
function open_threecylinder(; testnocon=false)

    # getting connection
    if testnocon
        dev = stepper_open(6, testnocon=true)
    else
        dev = stepper_open(6)
    end 

    # transformation functions cylinder 1
    θ0₁ = asin(-41.5/80.0)
    s2cx1(sx) = 80*(sin(θ0₁+ 2π*sx/2099) - sin(θ0₁))
    s2cy1(s) = 2π*12.5*s[2]/2009 - 80*(cos(θ0₁ + 2*π*s[1]/2099) - cos(θ0₁)) + NNyx1(s[1])
    c2sx1(x) = Int(round( (2099/(2π))*(asin(x/80+sin(θ0₁)) - θ0₁) ))    
    #c2sy1(c) = Int(round( (2009/(2π*12.5))*(c[2] + 80*(cos(asin(c[1]/80+sin(θ0₁)))-cos(θ0₁))) ))
    approx_c2sy1(c) = (2009/(2π*12.5))*(c[2] + 80*(cos(asin(c[1]/80+sin(θ0₂)))-cos(θ0₂)))
    function c2sy1(c)
        sx = c2sx1(c[1])
        minval(sy) = sum(abs2, c[2] - s2cy1([sx,sy]))
        r = optimize(v -> minval(first(v)), [approx_c2sy1(c)], LBFGS())
        return Int(round(r.minimizer[1]))
    end


    # transformation functions cylinder 2
    θ0₂ = asin(-59.0/80.0)
    s2cx2(sx) = 80*(sin(θ0₂+ 2π*sx/2176) - sin(θ0₂))
    s2cy2(s) = 2π*12.5*s[2]/1985 - 80*(cos(θ0₂ + 2*π*s[1]/2176) - cos(θ0₂)) + NNyx2(s[1])
    c2sx2(x) = Int(round( (2176/(2π))*(asin(x/80+sin(θ0₂)) - θ0₂) ))    
    approx_c2sy2(c) = (1985/(2π*12.5))*(c[2] + 80*(cos(asin(c[1]/80+sin(θ0₂)))-cos(θ0₂)))
    function c2sy2(c)
        sx = c2sx2(c[1])
        minval(sy) = sum(abs2, c[2] - s2cy2([sx,sy]))
        r = optimize(v -> minval(first(v)), [approx_c2sy2(c)], LBFGS())
        return Int(round(r.minimizer[1]))
    end

    # transformation functions cylinder 3
    θ0₃ = asin(-75.0/80.0)
    s2cx3(sx) = 80*(sin(θ0₃+ 2π*sx/2097) - sin(θ0₃))
    s2cy3(s) = 2π*12.5*s[2]/2041 - 80*(cos(θ0₃ + 2*π*s[1]/2097) - cos(θ0₃)) + NNyx3(s[1])
    c2sx3(x) = Int(round( (2097/(2π))*(asin(x/80+sin(θ0₃)) - θ0₃) ))    
    approx_c2sy3(c) = (2041/(2π*12.5))*(c[2] + 80*(cos(asin(c[1]/80+sin(θ0₃)))-cos(θ0₃)))
    function c2sy3(c)
        sx = c2sx3(c[1])
        minval(sy) = sum(abs2, c[2] - s2cy3([sx,sy]))
        r = optimize(v -> minval(first(v)), [approx_c2sy3(c)], LBFGS())
        return Int(round(r.minimizer[1]))
    end

    stepper_config!(dev, motorID=["x1","y1","x2","y2","x3","y3"], step2coord=[s2cx1,s2cy1,s2cx2,s2cy2,s2cx3,s2cy3], 
                    coord2step=[c2sx1,c2sy1,c2sx2,c2sy2,c2sx3,c2sy3], depend=[1,[1,2],3,[3,4],5,[5,6]])

    return dev
end



"""
    movecyl!(dev::StepperSystem, posnew::OneCylType; relat=false, strategy=1:2)

## Description
Function to move the cylinders system
## Arguments
- dev: stepper system from StepperControl library
- posnew: new CylType position
- strategy: stepper motor trigger order. If missing, then beststrategy() is called.
- relat: if true then the relative movement is performed. Default=false.
- method: method provided by the StepperControl library
"""
function movecyl!(dev::StepperSystem, posnew::CylType; strategy=missing, 
                    relat=false, method=StepperControl.manhattan_msg)

    coords = convert(Array, posnew)
    
    if typeof(posnew) == TwoCylType
        coords[3] += coords[1]
        coords[4] += coords[2]
    elseif typeof(posnew) == ThreeCylType
        coords[3] += coords[1]
        coords[4] += coords[2]
        coords[5] += coords[3]
        coords[6] += coords[4]
    end

    if ismissing(strategy)
        strategy = beststrategy(dev,posnew)
    end
    
    msg = stepper_move!(dev, coords, relat=relat, order=strategy, method=method)

    println(msg)
end
