#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado 
#                   and Prof. Sérgio Roberto Lopes
#       UFPR - 2024
#
#       This script works like a library that calculates a
#   recurrence plot from a time series and can use it to
#   calculate the probability of microstates and their 
#   entropy.
#
# ----------------------------------------------------------------------------------------- #
#       ** Required libraries:                                                              #
#   - Distances.jl                                                                          #
#   - StatsBase.jl                                                                          #
#   - Statistics.jl                                                                         #
#   - LinearAlgebra.jl                                                                      #
# ----------------------------------------------------------------------------------------- #
#       ** How do I use it?                                                                 #
#   - Include the file in your project using `include('../rp.jl')`.                         #
#   - Say that you want to use the RP module with `using .RP`.                              #
# ----------------------------------------------------------------------------------------- #
#       ** To do:                                                                           #
#   - Include the module in the Julia repository so that anyone can                         #
#   download it using Pkg. (I don't know how to do it >.<)                                  #
# ----------------------------------------------------------------------------------------- #
module RP                                                                                   #
# ----------------------------------------------------------------------------------------- #
using Distances                                                                             #
using StatsBase                                                                             #
using Statistics                                                                            #
using LinearAlgebra                                                                         #

using CairoMakie
# ----------------------------------------------------------------------------------------- #
export RecurrenceMatrix                                                                     #
"""
        Computes a recurrence matrix from a time series using a defined threshold ε and a
        recurrence function rrc (we use Θ(ε - |x - y|) by default).
"""
function RecurrenceMatrix(x, ε; rrc=std_recurrence)                                         #
    rp_size = size(x, 1)                                                                    #
    rp = zeros(Int, rp_size, rp_size)                                                       #
    #                                                                                       #
    for j = 1:rp_size                                                                       #
        for i = 1:j                                                                         #
            rp[i, j] = rrc(x[i, :], x[j, :], ε)                                             #
        end                                                                                 #
    end                                                                                     #
    #                                                                                       #
    return Symmetric(rp)                                                                    #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
export MicrostatesProbability                                                               #
"""
        It splits one recurrence structure into small parts like tiles and calculates the
        probability of finding each different tile in our structure. We call this 
        probability distribution the microstate probability.
"""
function MicrostatesProbability(rp, n::Int; samples=(0, 0),                                 #
    powvec::Vector=power_vector(n), no_use_samples=false)                                   #
    # ------------------------------------------------------------------------------------- #
    #           We use a dict here to save memory when we don't have all the microstates    #
    #   in our system. This works very well if we want to use n greater than 3 >.<          #
    tiles = Dict{Int,Float64}()                                                             #
    # ------------------------------------------------------------------------------------- #
    p = 0                                                                                   #
    add = 0                                                                                 #
    a_bin = 0                                                                               #
    count_samples = 0                                                                       #
    rp_size = size(rp)                                                                      #
    # ------------------------------------------------------------------------------------- #
    #           I will to verify here if the user want to use samples ...                   #
    index_x = missing                                                                       #
    index_y = missing                                                                       #
    if (no_use_samples)                                                                     #
        index_x = 1:n:rp_size[1]-(n-1)                                                      #
        index_y = 1:n:rp_size[2]-(n-1)                                                      #
    else                                                                                    #
        index_x = sample(1:rp_size[1]-(n-1), samples[1])                                    #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    for i in index_x                                                                        #
        if (!no_use_samples)                                                                #
            index_y = sample(1:rp_size[2]-(n-1), samples[2])                                #
        end                                                                                 #
        # --------------------------------------------------------------------------------- #
        for j in index_y                                                                    #
            add = 0                                                                         #
            p = 0                                                                           #
            #                                                                               #
            for x = 1:n                                                                     #
                for y = 1:n                                                                 #
                    a_bin = rp[i+x-1, j+y-1] == true ? 1 : 0                                #
                    add = add + a_bin * powvec[y+((x-1)*n)]                                 #
                end                                                                         #
            end                                                                             #
            #                                                                               #
            count_samples += 1                                                              #
            p = Int64(add) + 1                                                              #
            tiles[p] = get(tiles, p, 0) + 1                                                 #
        end                                                                                 #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    for index in keys(tiles)                                                                #
        tiles[index] /= (1.0 * count_samples)                                               #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    return tiles, count_samples                                                             #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
export Entropy                                                                              #
"""
        It calculates the entropy from the probabilities using Shannon's entropy.
"""
function Entropy(probs)                                                                     #
    S = 0                                                                                   #
    for index in keys(probs)                                                                #
        if (probs[index] > 0)                                                               #
            S += (-probs[index] * log(probs[index]))                                        #
        end                                                                                 #
    end                                                                                     #
    #                                                                                       #
    return S                                                                                #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
        To calculate the microstates we need a power vector to write them as
        binary numbers and convert them to decimal so that we can use it as an index to our dict. 
        So this function creates a vector where each element is a binary base 2^x where x
        goes from 0 to (n * n) - 1.
"""
function power_vector(n::Int)                                                               #
    powvec = zeros(Int64, (n * n))                                                          #
    for i in eachindex(powvec)                                                              #
        powvec[i] = Int64(2^(i - 1))                                                        #
    end                                                                                     #
    return powvec                                                                           #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
        Calculate the recurrence between two points in an n-dim space using: 
        r = Θ(ε - |x - y|)
"""
function std_recurrence(x, y, ε::Float64)                                                   #
    d = ε - euclidean(x, y)                                                                 #
    return Θ(d)                                                                             #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
        Calculate the recurrence between two points in an n-dim space using: 
        r = Θ(ε_max - |x - y|) * Θ(|x - y| - ε_min)
"""
function crr_recurrence(x, y, ε)                                                            #
    d = euclidean(x, y)                                                                     #
    return Θ(ε[2] - d) * Θ(d - ε[1])                                                        #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
        Calculates the threshold that maximises the recurrence entropy of the system.
        This function works to std recurrence: r = Θ(ε - |x - y|)
"""
function threshold_for_max_entroy_std_rrc(data, n; no_use_random_tiles=false,               #
    init_ε=0.00000001, ends_ε=0.49999999, macro_accuracy=10, micro_accurary=20)             #
    powvec = power_vector(n)                                                                #
    data_sz = size(data)                                                                    #
    # ------------------------------------------------------------------------------------- #
    max_s = 0.0                                                                             #
    max_threshold = ones(Float64, 3, data_sz[3]) * -1                                       #
    ε_macro_interval = range(init_ε, ends_ε, macro_accuracy)                                #
    # ------------------------------------------------------------------------------------- #
    #           First find the threshold that maximises the entropy in a largest interval.  #
    #   When we find it, we will run for a best accuracy in a smaller interval.             #
    for sample_index = 1:data_sz[3]                                                         #
        max_s = 0.0                                                                         #
        # --------------------------------------------------------------------------------- #
        for ε in ε_macro_interval                                                           #
            rp = RecurrenceMatrix(data[:, :, sample_index], ε)                              #
            probs, __tiles = MicrostatesProbability(rp, n; powvec=powvec,                   #
                samples=(floor(Int64, data_sz[1] / 4), floor(Int64, data_sz[1] / 4)),       #
                no_use_samples=no_use_random_tiles)                                         #
            s = Entropy(probs)                                                              #
            # ----------------------------------------------------------------------------- #
            if (s > max_s)                                                                  #
                max_threshold[1, sample_index] = max_threshold[2, sample_index]             #
                max_threshold[2, sample_index] = ε                                          #
                max_s = s                                                                   #
            else                                                                            #
                max_threshold[3, sample_index] = ε                                          #
                break                                                                       #
            end                                                                             #
        end                                                                                 #
        # --------------------------------------------------------------------------------- #
        max_s = 0.0                                                                         #
        # --------------------------------------------------------------------------------- #
        if (max_threshold[3, sample_index] < 0.0)                                           #
            throw("The maximum entropy wasn't found, try adjusting the thresholds or your data values.")
        end                                                                                 #
        # --------------------------------------------------------------------------------- #
        accuracy_interval = range(max_threshold[3, sample_index],                           #
            max_threshold[1, sample_index], 2 * micro_accurary)                             #
        for ε in accuracy_interval                                                          #
            rp = RecurrenceMatrix(data[:, :, sample_index], ε)                              #
            probs, __tiles = MicrostatesProbability(rp, n; powvec=powvec,                   #
                samples=(floor(Int64, data_sz[1] / 4), floor(Int64, data_sz[1] / 4)),       #
                no_use_samples=no_use_random_tiles)                                         #
            s = Entropy(probs)                                                              #
            # ----------------------------------------------------------------------------- #
            if (s > max_s)                                                                  #
                max_threshold[1, sample_index] = max_threshold[2, sample_index]             #
                max_threshold[2, sample_index] = ε                                          #
                max_s = s                                                                   #
            else                                                                            #
                max_threshold[3, sample_index] = ε                                          #
                break                                                                       #
            end                                                                             #
        end                                                                                 #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    return mean(max_threshold[2, :])                                                        #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
        Calculates the threshold that maximises the recurrence entropy of the system.
        This function works to std recurrence: r = Θ(ε_max - |x - y|) * Θ(|x - y| - ε_min)
"""
function threshold_for_max_entroy_crr_rrc(data, n; no_use_random_tiles=false,               #
    init_ε=0.00000001, ends_ε=0.49999999, macro_accuracy=10, micro_accurary=20)             #
    powvec = power_vector(n)                                                                #
    data_sz = size(data)                                                                    #
    # ------------------------------------------------------------------------------------- #
    max_threshold = zeros(Float64, data_sz[3], 2)                                           #
    ε_macro_interval = range(init_ε, ends_ε, macro_accuracy)                                #

    etr_space = zeros(Float64, macro_accuracy, macro_accuracy)
    etr_accur = zeros(Float64, micro_accurary, micro_accurary)

    # ------------------------------------------------------------------------------------- #
    for sample_index = 1:data_sz[3]                                                         #
        # --------------------------------------------------------------------------------- #
        for ε_max_index in eachindex(ε_macro_interval)                                      #
            for ε_min_index = 1:ε_max_index-1                                               #
                rp = RecurrenceMatrix(data[:, :, sample_index],                             #
                    (ε_macro_interval[ε_min_index], ε_macro_interval[ε_max_index]);         #
                    rrc=crr_recurrence)                                                     #
                probs, __tiles = MicrostatesProbability(rp, n; powvec=powvec,               #
                    samples=(floor(Int64, data_sz[1] / 4), floor(Int64, data_sz[1] / 4)),   #
                    no_use_samples=no_use_random_tiles)                                     #
                s = Entropy(probs)                                                          #
                # ------------------------------------------------------------------------- #
                etr_space[ε_min_index, ε_max_index] = s                                     #
            end                                                                             #
        end                                                                                 #
        # --------------------------------------------------------------------------------- #
        #       Pega os arredores do ponto de máximo.
        mx = findmax(etr_space)
        # --------------------------------------------------------------------------------- #
        #       Faz o intervalo de precissão
        ε_min_accuracy_interval = range(ε_macro_interval[mx[2][1]-1],
            ε_macro_interval[mx[2][1]+1], micro_accurary)
        ε_max_accuracy_interval = range(ε_macro_interval[mx[2][2]-1],
            ε_macro_interval[mx[2][2]+1], micro_accurary)
        # --------------------------------------------------------------------------------- #
        for ε_max_index in eachindex(ε_max_accuracy_interval)                               #
            for ε_min_index = eachindex(ε_min_accuracy_interval)                            #
                rp = RecurrenceMatrix(data[:, :, sample_index],                             #
                    (ε_min_accuracy_interval[ε_min_index],                                  #
                        ε_max_accuracy_interval[ε_max_index]); rrc=crr_recurrence)          #
                probs, __tiles = MicrostatesProbability(rp, n; powvec=powvec,               #
                    samples=(floor(Int64, data_sz[1] / 4), floor(Int64, data_sz[1] / 4)),   #
                    no_use_samples=no_use_random_tiles)                                     #
                s = Entropy(probs)                                                          #
                # ------------------------------------------------------------------------- #
                etr_accur[ε_min_index, ε_max_index] = s                                     #
            end                                                                             #
        end                                                                                 #
        # --------------------------------------------------------------------------------- #
        #       Pega o valor máximo da precissão                                            #
        mx_accr = findmax(etr_accur)                                                        #
        max_threshold[sample_index, 1] = ε_min_accuracy_interval[mx_accr[2][1]]
        max_threshold[sample_index, 2] = ε_max_accuracy_interval[mx_accr[2][2]]
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    return (mean(max_threshold[:, 1]), mean(max_threshold[:, 2]))
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
    Function heaviside.
"""
function Θ(x::Float64)                                                                      #
    if (x >= 0)                                                                             #
        return 1                                                                            #
    end                                                                                     #
    return 0                                                                                #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
end                                                                                         #