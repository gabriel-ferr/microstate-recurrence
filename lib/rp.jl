#
#       By Gabriel Ferreira
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
export RecurrenceMatrix
# ------------------------------------------------------------------
"""
        RecurrenceMatrix(x[, dims], ε::Any (1); rrc = StdRrc)
    
    (1) - type of variable associated with the threshold of the 
    recurrence function (rrc).

    Use a user-defined recurrence function to compute and return a
    recurrence matrix from a time series.
"""
function RecurrenceMatrix(x::Array{Float64,2}, ε::Any; rrc='a')
    rp_size = size(x, 1)
    rp = zeros(Int, rp_size, rp_size)

    for j = 1:rp_size
        for i = 1:j
            rp[i, j] = rrc(x[i, :], x[j, :], ε)
        end
    end
end
# ------------------------------------------------------------------
end