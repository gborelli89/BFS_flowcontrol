# -----------------------------------------------------------------------
# Data structure
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Cylinder position data structure
# -----------------------------------------------------------------------

# The two first values refer to a global coordinate system 

# All the other values are local with the origin on the upwind cylinder center

# Global origin at the step edge

# These were made mutable in order to be changed without creating 
# new variables, specially to get the current position.

# -----------------------------------------------------------------------
# Required packages
using AutoHashEquals, DataFrames

# 1 cylinder position
@auto_hash_equals mutable struct OneCylType
    x::Float64 # cylinder center - x coordinate (global)
    y::Float64 # cylinder center - y coordinate (global)
end



# 2 cylinder position
@auto_hash_equals mutable struct TwoCylType
    x::Float64 # first cylinder center - x coordinate (global)
    y::Float64 # first cylinder center - y coordinate (global)
    dx::Float64 # delta x to the second cylinder (related to the first)
    dy::Float64 # delta y to the second cylinder (related to the first)
end


# 3 cylinder position
@auto_hash_equals mutable struct ThreeCylType
    x::Float64 # first cylinder center - x coordinate (global)
    y::Float64 # first cylinder center - y coordinate (global)
    dx1::Float64 # delta x to the second cylinder (related to the first)
    dy1::Float64 # delta y to the second cylinder (related to the first)
    dx2::Float64 # delta x to the third cylinder (related to the second)
    dy2::Float64 # delta y to the third cylinder (related to the second)
end

# Union
CylType = Union{OneCylType,TwoCylType,ThreeCylType}
ArrayCylType = Union{Array{OneCylType,1},Array{TwoCylType,1},Array{ThreeCylType,1}}

# -----------------------------------------------------------------------
# Data structure for the stepper motor
# -----------------------------------------------------------------------
# No need for mutable data here.
# -----------------------------------------------------------------------

# Step data and strategy of mov. to pass to the stepper motor
@auto_hash_equals struct StepperMotionType
    step::Array{Int64,1}
    strategy::Array{Int64,1}
end

# -----------------------------------------------------------------------
# New method for convert
# -----------------------------------------------------------------------
# Make convert work for all the structures created
# -----------------------------------------------------------------------

# For one cylinder
function Base.convert(::Type{Array}, x::OneCylType)
    return [x.x, x.y]
end

function Base.convert(::Type{DataFrame}, x::OneCylType)
    return DataFrame(x=x.x, y=x.y)
end

function Base.convert(::Type{DataFrame}, x::Array{OneCylType,1})
    v = DataFrame(x=Float64[],y=Float64[])
    for i in x
        append!(v, convert(DataFrame, i))
    end
    return v
end

function Base.convert(::Type{OneCylType}, x::AbstractVector)
    if length(x) != 2
        throw(DimensionMismatch("vector length should be equal 2"))
    end
    return OneCylType(x[1],x[2])
end

function Base.convert(::Type{OneCylType}, x::DataFrame)
    if nrow(x) == 1
        return OneCylType(x.x, x.y)
    else
        e = eachrow(x)
        return [OneCylType(i.x, i.y) for i in e]
    end
end

# For two cylinder
function Base.convert(::Type{Array}, x::TwoCylType)
    return [x.x, x.y, x.dx, x.dy]
end

function Base.convert(::Type{DataFrame}, x::TwoCylType)
    return DataFrame(x=x.x, y=x.y, dx=x.dx, dy=x.dy)
end

function Base.convert(::Type{DataFrame}, x::Array{TwoCylType,1})
    v = DataFrame(x=Float64[],y=Float64[],dx=Float64[],dy=Float64[])
    for i in x
        append!(v, convert(DataFrame, i))
    end
    return v
end

function Base.convert(::Type{TwoCylType}, x::AbstractVector)
    if length(x) != 4
        throw(DimensionMismatch("vector length should be equal 4"))
    end
    return TwoCylType(x[1],x[2],x[3],x[4])
end

function Base.convert(::Type{TwoCylType}, x::DataFrame)
    if nrow(x) == 1
        return TwoCylType(x.x, x.y, x.dx, x.dy)
    else
        e = eachrow(x)
        return [TwoCylType(i.x, i.y, i.dx, i.dy) for i in e]
    end
end

# For three cylinder
function Base.convert(::Type{Array}, x::ThreeCylType)
    return [x.x, x.y, x.dx1, x.dy1, x.dx2, x.dy2]
end

function Base.convert(::Type{DataFrame}, x::ThreeCylType)
    return DataFrame(x=x.x, y=x.y, dx1=x.dx1, dy1=x.dy1, dx2=x.dx2, dy2=x.dy2)
end

function Base.convert(::Type{DataFrame}, x::Array{ThreeCylType,1})
    v = DataFrame(x=Float64[],y=Float64[],dx1=Float64[],dy1=Float64[],dx2=Float64[],dy2=Float64[])
    for i in x
        append!(v, convert(DataFrame, i))
    end
    return v
end

function Base.convert(::Type{ThreeCylType}, x::AbstractVector)
    if length(x) != 6
        throw(DimensionMismatch("vector length should be equal 6"))
    end
    return ThreeCylType(x[1],x[2],x[3],x[4],x[5],x[6])
end

function Base.convert(::Type{ThreeCylType}, x::DataFrame)
    if nrow(x) == 1
        return ThreeCylType(x.x, x.y, x.dx1, x.dy1, x.dx2, x.dy2)
    else
        e = eachrow(x)
        return [ThreeCylType(i.x, i.y, i.dx1, i.dy1, i.dx2, i.dy2) for i in e]
    end
end

# -----------------------------------------------------------------------
# New method for getindex 
# -----------------------------------------------------------------------
# Make getindex work for all the structures created
# -----------------------------------------------------------------------
function Base.getindex(cyl::CylType, i::Integer)
    a = convert(Array, cyl)
    return a[i]
end

function Base.getindex(cyl::CylType, r::AbstractRange)
    a = convert(Array, cyl)
    return a[r]
end


# For StepperMotionType
function Base.getindex(sm::StepperMotionType, i::Integer)
    if i==1
        sm.step
    elseif i==2
        sm.strategy
    else
        throw(DimensionMismatch("Dimensions not matching"))
    end
end

# -----------------------------------------------------------------------
# New method for setindex! 
# -----------------------------------------------------------------------
# Make getindex work for all the structures created
# -----------------------------------------------------------------------

# For OneCylType
function Base.setindex!(cyl::OneCylType, v::AbstractFloat, i::Integer)
    if i == 1
        cyl.x = v
    elseif i == 2
        cyl.y = v
    else
        throw(BoundsError)
    end
end

# For TwoCylType
function Base.setindex!(cyl::TwoCylType, v::AbstractFloat, i::Integer)
    if i==1
        cyl.x = v
    elseif i==2
        cyl.y = v
    elseif i==3
        cyl.dx = v
    elseif i==4
        cyl.dy = v
    else
        throw(BoundsError)
    end
end

# For ThreeCylType
function Base.setindex!(cyl::ThreeCylType, v::AbstractFloat, i::Integer)
    if i==1
        cyl.x = v
    elseif i==2
        cyl.y = v
    elseif i==3
        cyl.dx1 = v
    elseif i==4
        cyl.dy1 = v
    elseif i == 5
        cyl.dx2 = v
    elseif i == 6
        cyl.dy2 = v
    else
        throw(BoundsError)
    end
end

# -----------------------------------------------------------------------
# New method for length
# -----------------------------------------------------------------------
# Make length work for all the structures created
# -----------------------------------------------------------------------
Base.length(cyl::CylType) = convert(Array, cyl)


# -----------------------------------------------------------------------
# New method for length
# -----------------------------------------------------------------------
# Make length work for all the structures created
# -----------------------------------------------------------------------
Base.lastindex(cyl::CylType) = length(cyl)


Base.length(sm::StepperMotionType) = 2

