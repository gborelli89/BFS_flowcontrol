# -----------------------------------------------------------------------
# Movement optimization
# -----------------------------------------------------------------------
# OBS.:
# These algorithms were not during the tests! 
# It was obseved little impact on the experiments runtime.
# If necessary the user can adapt and use it to determine moving order
# -----------------------------------------------------------------------
# Compute distances 
# -----------------------------------------------------------------------
using Distances

"""
    distcylinder(posOld::CylType, posNew::CylType; metric="steps")

## Description
Compute distance between two positions.
## Arguments
- posOld: previous position (from)
- posNew: new position (to)
- metric: metric to compute distance: "euclidean", "manhattan", "steps"
# Examples
```jldoctest
julia> distcylinder(OneCylType(0,0), OneCylType(10,0))
57

julia> distcylinder(OneCylType(0,0), OneCylType(10,10), metric="manhattan")
20.0

julia> distcylinder(ThreeCylType(0,0,5,0,5,0), ThreeCylType(10,5,5,0,5,1))
696
```
"""
function distcylinder(posOld::OneCylType, posNew::OneCylType; metric="steps")
    if metric == "euclidean"
        r = euclidean([posNew.x,posNew.y], [posOld.x,posOld.y])
    elseif metric == "manhattan"
        r = abs(posNew.x-posOld.x) + abs(posNew.y-posOld.y)
    elseif metric == "steps"
        steps = coord2steps(posOld, posNew)[1].step
        r = sum(abs.(steps))
    else
        error("Metric is not valid")
    end
    return r
end

function distcylinder(posOld::TwoCylType, posNew::TwoCylType; metric="steps")
    if metric == "euclidean"
        r1 = euclidean([posNew.x,posNew.y], [posOld.x,posOld.y])
        r2 = euclidean([posNew.x+posNew.dx, posNew.y+posNew.dy], 
                        [posOld.x+posOld.dx, posOld.y+posOld.dy])
        r = r1 + r2
    elseif metric == "manhattan"
        r1 = abs(posNew.x-posOld.x) + abs(posNew.y-posOld.y)
        r2 = r1 + abs(posNew.dx-posOld.dx) + abs(posNew.dy-posOld.dy)
        r = r1 + r2
    elseif metric == "steps"
        steps = coord2steps(posOld, posNew)[1].step
        r = sum(abs.(steps))
    else
        error("Metric is not valid")
    end
    return r
end

function distcylinder(posOld::ThreeCylType, posNew::ThreeCylType; metric="steps")
    if metric == "euclidean"
        r1 = euclidean([posNew.x,posNew.y], [posOld.x,posOld.y])
        r2 = euclidean([posNew.x+posNew.dx1,posNew.y+posNew.dy1], 
                        [posOld.x+posOld.dx1,posOld.y+posOld.dy1])
        r3 = euclidean([posNew.x+posNew.dx1+posNew.dx2,posNew.y+posNew.dy1+posNew.dy2], 
                        [posOld.x+posOld.dx1+posOld.dx2,posOld.y+posOld.dy1+posOld.dy2])
        r = r1 + r2 + r3
    elseif metric == "manhattan"
        r1 = abs(posNew.x-posOld.x) + abs(posNew.y-posOld.y)
        r2 = r1 + abs(posNew.dx1-posOld.dx1) + abs(posNew.dy1-posOld.dy1)
        r3 = r2 + abs(posNew.dx2-posOld.dx2) + abs(posNew.dy2-posOld.dy2)
        r = r1 + r2 + r3
    elseif metric == "steps"
        steps = coord2steps(posOld, posNew)[1].step
        r = sum(abs.(steps))
    end
    return r
end

"""
    greedypath(actualPos::CylType, points::ArrayCylType)

## Description
Greedy algorithm to determine the near optimum path (nearest neighbors).
## Arguments
- actualPos: the initial position of the robot
- points: positions the robot must visit
"""
function greedypath(actualPos::CylType, points::ArrayCylType)

    if typeof(actualPos) != typeof(points[1])
        error("actualPos and points don't have elements with the same type")
    end

    path = [] # initialize path

    n = length(points)
    notVisited = repeat([true],n)
    id = collect(1:n)

    while sum(notVisited) > 0
        minDist = argmin([distcylinder(actualPos, i) for i in points[notVisited]])
        idtemp = id[notVisited][minDist]
        path = push!(path, idtemp)
        notVisited[idtemp] = false
        actualPos = points[idtemp]
    end

    return path
end 

"""
    benchoptimizer(n=50; nc=3, optimFun=greedypath)

## Description
Benchmark tool for the path optimizer (compared to random ordered).
## Arguments
- n: number of points
- nc: number of cylinders
- fun: optimization function (only greeypath allowed yet)
"""
function benchoptimizer(n=50; nc=3, fun=greedypath)
    
    # One cylinder
    initValue = randomInd(nc)
    p = initPop(n,nc=nc)
    s = distcylinder(initValue, p[1])

    id_optim = fun(initValue,p)
    s_optim = distcylinder(initValue, p[id_optim[1]])
    
    for i in 2:n
        k = i-1
        s += distcylinder(p[k], p[i])
        s_optim += distcylinder(p[id_optim[k]], p[id_optim[i]])
    end

    println("Random order: ", s, " steps")
    println("Optimized: ", s_optim, " steps")
    println("Ratio: ", s_optim/s)

    return s_optim/s
end
