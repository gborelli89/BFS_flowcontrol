# Only used for simulations of one cylinder in ga_binary.jl

using Random
using DataFrames
using CSV

# Generate position Array for one cylinder
# x,y: coordinate positions
function generatePosArray(x, y)
    p = Array{OneCylType,1}()
    for i in x
        for j in y
            push!(p, OneCylType(i,j))
        end
    end
    return p
end

# Get test data for one cylinder
function simDataOne(dname; xoffset=0, yoffset=0)
    x = DataFrame(CSV.read(dname*"/x.csv", header=false))
    x = collect(x[:,1]) .+ xoffset
    y = DataFrame(CSV.read(dname*"/y.csv", header=false))
    y = collect(y[:,1]) .+ yoffset
    data = DataFrame(CSV.read(dname*"/pressureData.csv", header=1))

    dict = Dict{OneCylType, Float64}()

    for i in 1:length(x)
        for k in 1:length(y)
            dict[OneCylType(x[i],y[k])] = data[i,k]
        end
    end

    return dict

end


# Generate randomic data for simulation (not working well)
# nc: number of cylinders
function simDataRandom(nc)
    if nc == 1 # OneCylType for size=[6,5]
        x = 0:1:63
        y = -15:1:16
        p = generatePosArray(x,y)
        n = length(p)

        x1 = 10:1:60
        y1 = -5:1:8
        p1 = generatePosArray(x1,y1)
        n1 = length(p1)

        x2 = 20:1:25
        y2 = 3:1:5
        p2 = generatePosArray(x2,y2)
        n2 = length(p2)

        d = Dict{OneCylType, Float64}()
        [d[p[k]] = 10 + rand(Random.seed!(k)) for k in 1:n]
        [d[p1[k]] = d[p1[k]] + 2 + rand(Random.seed!(k)) for k in 1:n1]
        [d[p2[k]] = d[p2[k]] + 2 + 2*rand(Random.seed!(k)) for k in 1:n2]
    else
        error("Dimension mismatch!")
    end
    return d
end

