#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#       UFPR - 2024
#
#       This script works like a library that calculates a
#   recurrence plot from a time series and can use it to
#   calculate the probability of microstates and their 
#   entropy.
# ------------------------------------------------------------------
#       ** Required libraries:
#   - Distances.jl
#   - StatsBase.jl
#   - Statistics.jl
#   - CairoMakie.jl
#   - LinearAlgebra.jl
# ------------------------------------------------------------------
#       ** How do I use it?
#   - Include the file in your project using `include('../rp.jl')`.
#   - Say that you want to use the RP module with `using .RP`.
# ------------------------------------------------------------------
#       ** To do:
#   - Include the module in the Julia repository so that anyone can
#   download it using Pkg. (I don't know how to do it >.<)
# ------------------------------------------------------------------
module RP
# ------------------------------------------------------------------
#   Libraries:
using Distances
using StatsBase
using Statistics
using CairoMakie
using LinearAlgebra
# ------------------------------------------------------------------
#   Exports (functions that can be used in an external environment.)
export Θ
export StdRrc
export CrrRrc
export Entropy
export RecurrencePlot
export GetMicrostates
export GetPowerVector
export RecurrenceMatrix
# ------------------------------------------------------------------
"""
    RecurrenceMatrix(x[, dims], ε::Any (1); rrc = StdRrc)
    
(1) - type of variable associated with the threshold of the 
recurrence function (rrc).

Use a user-defined recurrence function to compute and return a
recurrence matrix from a time series.
"""
function RecurrenceMatrix(x, ε::Any; rrc=StdRrc)
    rp_size = size(x, 1)
    rp = zeros(Int, rp_size, rp_size)

    for j = 1:rp_size
        for i = 1:j
            rp[i, j] = rrc(x[i, :], x[j, :], ε)
        end
    end

    return Symmetric(rp)
end
# ------------------------------------------------------------------
"""
        RecurrencePlot(RecurrenceMatrix; xlabel="Time", ylabel="Time")

    Uses the CairoMakie to plot a RecurrencePlot from a RecurrenceMatrix.
"""
function RecurrencePlot(rp; xlabel="Time", ylabel="Time")
    fig = Figure()
    ax = Axis(fig[1, 1])
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    heatmap!(ax, rp, colormap=:binary)

    return fig
end
# ------------------------------------------------------------------
"""
    GetMicrostates(RecurrenceMatrix, n::Int; powvec::Vector = GetPowerVector(n::Int) (1))

(1) If you are going to run this function in a loop, it is recommended
that you give it a pre-instantiated power vector (powvec), use
function GetPowerVector(n::Int) to do it.

This function calculates the microstates of form nxn from a recurrence
matrix given to it and returns the probability of each microstate 
being located.

"""
function GetMicrostates(rp, n::Int; powvec::Vector=GetPowerVector(n))
    #   It uses a dict to save memory when we don't have all the
    #   microstates in our system.
    R = Dict{Int,Float64}()
    # --------------------------------------------------------------
    p = 0
    add = 0
    a_bin = 0
    count_samples = 0
    rp_size = size(rp)
    # --------------------------------------------------------------
    for j = 1:n:rp_size[2]-(n-1)
        for i = 1:n:rp_size[1]-(n-1)
            add = 0
            p = 0
            for x = 1:n
                for y = 1:n
                    a_bin = rp[i+x-1, j+y-1] == true ? 1 : 0
                    add = add + a_bin * powvec[y+((x-1)*n)]
                end
            end
            count_samples += 1
            p = Int64(add) + 1
            R[p] = get(R, p, 0) + 1
        end
    end
    # --------------------------------------------------------------
    for r in keys(R)
        R[r] /= (1.0 * count_samples)
    end
    # --------------------------------------------------------------
    return R
end
# ------------------------------------------------------------------
"""
    GetMicrostates(RecurrenceMatrix, n::Int, samples::Tuple(n_rows::Int, n_samples::Int); powvec = GetPowerVector(n::Int) (1))

(1) If you are going to run this function in a loop, it is recommended
that you give it a pre-instantiated power vector (powvec), use
function GetPowerVector(n::Int) to do it.

**  I haven't tested this function yet >.<

This function calculates the microstates of form nxn from a recurrence
matrix given to it and returns the probability of each microstate 
being located.

Sampling: it uses a sampling system that gets a unique collection of
indeces for each row and uses it to calculate the microstate 
probabilities. Note that n_rows must be less than size(p, 1) - (n - 1)
and the same thing occur with n_samples but with size(p, 2) - (n - 1).

"""
function GetMicrostates(rp, n::Int, samples; powvec=GetPowerVector(n))
    if (samples[1] > size(rp, 1) - (n - 1))
        throw("The number of n_rows exceeds the maximum.")
    end
    if (samples[2] > size(rp, 2) - (n - 1))
        throw("The number of n_samples exceeds the maximum.")
    end
    # --------------------------------------------------------------
    #   It uses a dict to save memory when we don't have all the
    #   microstates in our system.
    R = Dict{Int,Float64}()
    # --------------------------------------------------------------
    p = 0
    add = 0
    a_bin = 0
    rp_size = size(rp)
    index_x = sample(1:rp_size[1]-(n-1), samples[1])
    index_y = sample(1:rp_size[2]-(n-1), samples[2])
    # --------------------------------------------------------------
    for i in index_x
        #   For each row, take one sample.
        index_y = sample(1:rp_size[2]-(n-1), sample[2])
        for j in index_y
            add = 0
            p = 0
            for x = 1:n
                for y = 1:n
                    a_bin = rp[i+x-1, j+y-1] == true ? 1 : 0
                    add = add + a_bin * powvec[y+((x-1)*n)]
                end
            end
            p = Int64(add) + 1
            R[p] = get(R, p, 0) + 1
        end
    end
    # --------------------------------------------------------------
    for r in keys(R)
        R[r] /= (1.0 * (sample[1] * sample[2]))
    end
    # --------------------------------------------------------------
    return R
end
# ------------------------------------------------------------------
"""
    Entropy(probs::Dict{Int, Float64})

It calculates the entropy from the microstate probabilities using
Shannon's entropy.
"""
function Entropy(probs::Dict{Int,Float64})
    S = 0
    for r in keys(probs)
        if (probs[r] > 0)
            S += (-probs[r] * (log(probs[r])))
        end
    end
    return S
end
# ------------------------------------------------------------------
"""
    GetPowerVector(n::Int)

To calculate the microstates we need a power vector to write them as
binary numbers and use as an index to our dict. So this function
creates a vector where each element is a binary base 2^x where x
goes from 0 to (n * n) - 1.
"""
function GetPowerVector(n::Int)
    powvec = zeros(Int64, (n * n))
    for i in eachindex(powvec)
        powvec[i] = Int64(2^(i - 1))
    end
    return powvec
end
# ------------------------------------------------------------------
"""
    StdRrc(x::Vector{Float64}, y::Vector{Float64}, ε::Float64)

Calculate the recurrence between two points in an n-dim space 
using: r = Θ(ε - |x - y|)
"""
function StdRrc(x, y, ε::Float64)
    d = ε - euclidean(x, y)
    return Θ(d)
end
# ------------------------------------------------------------------
"""
    CrrRrc(x::Vector{Float64}, y::Vector{Float64}, ε::Tuple(ε_min::Float64, ε_max::Float64))

Calculate the recurrence between two points in an n-dim space
using: r = Θ(ε_max - |x - y|) * Θ(|x - y| - ε_min)
"""
function CrrRrc(x, y, ε)
    d = euclidean(x, y)
    return Θ(ε[2] - d) * Θ(d - ε[1])
end
# ------------------------------------------------------------------
"""
    Function heaviside.
"""
function Θ(x::Float64)
    if (x >= 0)
        return 1
    end
    return 0
end
# ------------------------------------------------------------------
"""
    MaximizeStdEntropy(x::Array{Float64, 3}, n::Int; samples::Int = 1000, ε_min::Float64 = 0, ε_max::Float64 = 1)

Compute the maximum entropy for a time serie when we use the Standard Recurrence.
"""
function MaximizeStdEntropy(x, n::Int; samples=1000, ε_min=0, ε_max=1)
    ε = range(ε_min, ε_max, samples)
    smat = zeros(Float64, samples)
    x_size = size(x)
    powvec = GetPowerVector(n)
    ε_res = zeros(Float64, x_size[3])
    # --------------------------------------------------------------
    for s in x_size[3]
        smat = zeros(Float64, samples)
        Threads.@threads for i in eachindex(ε)
            rp = RecurrenceMatrix(x[:, :, s], ε[i]; rrc=StdRrc)
            smat[i] = Entropy(GetMicrostates(rp, n; powvec=powvec))
        end
        ε_res[s] = ε[findmax(smat)[2]]
    end
    # --------------------------------------------------------------
    return mean(ε_res)
end
# ------------------------------------------------------------------
"""
    MaximizeStdEntropy(x::Array{Float64, 3}, n::Int; samples::Int = 1000, ε_min::Float64 = 0, ε_max::Float64 = 1)

Compute the maximum entropy for a time serie when we use the Corridor Recurrence.
"""
function MaximizeCrrEntropy(x, n::Int; samples=1000, ε_min=0, ε_max=1)
    ε = range(ε_min, ε_max, samples)
    smat = zeros(Float64, samples, samples)
    powvec = GetPowVector(n)
    x_size = size(x)
    ε_res = zeros(Float64, x_size[3], 2)
    # --------------------------------------------------------------
    for s in x_size[3]
        smat = zeros(Float64, samples, samples)
        Threads.@threads for i in eachindex(ε)
            for j = 1:(i-1)
                rp = RecurrenceMatrix(x[:, :, s], (ε[j], ε[i]); rrc=CrrRrc)
                smat[j, i] = Entropy(GetMicrostates(rp, n; powvec=powvec))
            end
        end
        mx = findmax(smat)[2]
        ε_res[s, 1] = ε[mx[1]]
        ε_res[s, 2] = ε[mx[2]]
    end
    # --------------------------------------------------------------
    return (mean(ε_res[:, 1]), mean(ε_res[:, 2]))
end
# ------------------------------------------------------------------
end