# Functions to perform stepper motor calibration
using LsqFit, JLD, StepperControl

"""
    create_xmodel(R, xorig)

## Description
Create a model for x displacements.
## Arguments
- R: bell crank length
- xorig: horizontal distance from the center of rotation to the origin
## Output
A function to convert steps, in x direction, into displacements with the following arguments:
- steps: number of steps (integer)
- p: parameters -> steps per revolution (to be find)
"""
function create_xmodel(R=80.0, xorig=-41.5)
  
    θ0 = asin(xorig/R)
    @. model(steps, p) = R*(sin(θ0 + 2π*steps/p[1]) - sin(θ0))
  
    return model
end


"""
    create_yxmodel(R, xorig)

## Description
Create a model for y displacements due to rotation of stepper motor responsible for x displacements.
## Arguments
- R: bell crank length
- xorig: horizontal distance from the center of rotation to the origin
## Output
A function to convert steps, in x direction, into displacements with the following arguments:
- steps: number of steps (integer)
- p: parameters -> steps per revolution (to be find)
"""
function create_yxmodel(R=80.0, xorig=-41.5)
  
    θ0 = asin(xorig/R)
    @. model(steps, p) = - R*(cos(θ0 + 2π*steps/p[1]) - cos(θ0))
  
    return model
end

"""
    create_ymodel(r)

## Description
Create a model for x displacements.
## Arguments
- r: gear radius
## Output
A function to convert steps, in y direction, into displacements with the following arguments:
- steps: number of steps (integer)
- p: parameters -> steps per revolution (to be find)
"""
function create_ymodel(r=12.5)
  
    @. model(steps, p) = 2π*steps*r/p[1]
  
    return model
end

"""
    iterative_calibr(dev::StepperSystem, steps::AbstractArray)

## Description
Iterative calibration process.
## Arguments
- dev: stepper motor device (see StepperControl)
- steps: array with the desired steps. 
## Obs
The step values must be given relative to the origin (absolute displacements).
"""
function iterative_calibr(dev::StepperSystem, steps::AbstractArray)

    measurement = Float64[]
    println("Press enter to start calibration.\n")
    readline()

    for i in steps
        msg = stepper_move!(dev, i, relat=false)
        println(msg)
        println("Enter displacement for ", i, " steps: ")
        push!(measurement, parse(Float64, readline()))
    end

    println("\nEnd of calibration!")
    return measurement
end


"""
    calibr1DOF(sname::AbstractString="calibr.jld"; id, model, p0, steps::AbstractArray; measurement=missing)

## Description
Calibration function.
## Arguments
- sname: filename to be saved
- id: stepper ID (see StepperControl)
- model: model to be fit
- p0: initial paramenters (according to the model, elements must be of float type)
- steps: array with the desired steps (absolute displacement will be applied)
- measurement: measured data for each step. If is missing, then a iterative calibration is called
"""
function calibr1DOF(sname::AbstractString="calibr.jld"; id, model, p0, steps::AbstractArray, measurement=missing)

    n = length(steps)
    steps = Int.(round.(steps))

    if ismissing(measurement)
        r = stepper_open(1)
        stepper_config!(r, motorID=[id]) 
        println("Open the serial monitor and press enter to start.")
        readline()   
        measurement = iterative_calibr(r, steps)
    end

    fit = curve_fit(model, steps, measurement, p0)
    save(sname, "fit", fit, "steps", steps, "measurement", measurement)

    return measurement, fit
end
