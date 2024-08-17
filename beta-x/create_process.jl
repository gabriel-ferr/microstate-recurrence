#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#       UFPR - 2024
#
#       Script to create a new process of execution before 
#   run the code.
# ------------------------------------------------------------------
include("../lib/rp.jl")
include("../lib/run_manager.jl")
# ------------------------------------------------------------------
using .RP
using .RManager
using Random
# ------------------------------------------------------------------
rng = MersenneTwister()
Random.seed!()
# ------------------------------------------------------------------
#       Calculate the time series with a fixed size after a 
#   transient.
function GenerateData(size, samples, β)
    data = zeros(Float64, (size, 1, samples, length(β)))
    x_init = rand(Float64, samples)
    transient = round(Int64, (10 * size))
    #---------------------------------------------------------------  
    for b in eachindex(β)
        for s in eachindex(x_init)
            x_before = x_init[s]
            for j = 1:(size+transient)
                x_after = x_before * β[b]
                while (x_after > 1.0)
                    x_after = x_after - 1.0
                end
                if (j > transient)
                    data[j-transient, 1, s, b] = x_after
                end
                x_before = x_after
            end
        end
    end
    #--------------------------------------------------------------- 
    return data
end
# ------------------------------------------------------------------
#   Size of the Beta X time series.
const β_size = 1000
#
#   Numb of samples used.
const β_samples = 5
#
#   Values of Beta X that we will use.
const β = [2.89, 3.49, 3.99, 4.59, 5.09, 6.19, 6.89, 7.59]
# ------------------------------------------------------------------
#       Values initial to x that we will use to calculate the 
#   maximum entropy.
etr_data = GenerateData(β_size, β_samples, β)
RManager.CreateTmpFile("etr_data", etr_data, ARGS[1])
