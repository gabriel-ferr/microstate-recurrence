#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado 
#                   and Prof. Sérgio Roberto Lopes
#       UFPR - 2024
#
#           This script creates a new process for the beta-x project. 
#   Basically, it defines a set of values for beta and an initial data 
#   set, and saves this data to a json file that can be restored if the
#   process needs to be stopped. This file also stores the progress of 
#   the run.
# ----------------------------------------------------------------------------------------- #
#       ** Required libraries:                                                              #
#   - Local: RManager.jl                                                                    #
#   - Random.jl                                                                             #
# ----------------------------------------------------------------------------------------- #
if !(@isdefined RManager)                                                                   #
    include("../lib/run_manager.jl")                                                        #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
using .RManager                                                                             #
using Random                                                                                #
# ----------------------------------------------------------------------------------------- #
rng = MersenneTwister()                                                                     #
Random.seed!()                                                                              #
# ----------------------------------------------------------------------------------------- #
#       ** Settings:    (Edit here to configure the generated process =D)                   #
#   - Length of the time series used in the operations.                                     #
const series_size = 1000                                                                    #
#                                                                                           #
#   - Number of samples used to calculate the maximum entropy mean.                         #
const max_entropy_samples = 5                                                               #
#                                                                                           #
#   - Number of samples used to train the network.                                          #
const net_train_samples = 2000                                                              #
#   - The values used to train the network must be random?                                  #
const net_train_ramdom = false                                                              #
#                                                                                           #
#   - Number of samples used to test the network.                                           #
#   These values will be generated to be different from the values used for training.       #
const net_test_samples = 1000                                                               #
#                                                                                           #
#   - Values of Beta X used =3                                                              #
const β = [2.89, 3.49, 3.99, 4.59, 5.09, 6.19, 6.89, 7.59]                                  #
#                                                                                           #
# ----------------------------------------------------------------------------------------- #
#           - Generates a set of beta x values                                              #
function generate_betax(size, samples; no_use_random=false, check_content=[])               #
    result = zeros(Float64, (size, 1, samples, length(β)))                                  #
    x_init = rand(Float64, samples)                                                         #
    transient = round(Int64, (10 * size))                                                   #
    # ------------------------------------------------------------------------------------- #
    if (length(check_content) > 0)                                                          #
        if (no_use_random)                                                                  #
            no_use_random = false                                                           #
            println("WARING: Using a checklist to generate data means that you can't use values that aren't random.")
        end                                                                                 #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    #       If you do not want the script to use random values, get them from a series      #
    #   between 0 and 1.                                                                    #
    if (no_use_random)                                                                      #
        x_init = range(0, 1, samples)                                                       #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    #       Guarantees that the x init set doesn't have any values listed in a check list.  #
    if (length(check_content) > 0)                                                          #
        for index in eachindex(x_init)                                                      #
            while (x_init[index] in check_content)                                          #
                new_value = rand(Float64, 1)                                                #
                while (new_value in x_init)                                                 #
                    new_value = rand(Float64, 1)                                            #
                end                                                                         #
                x_init[index] = new_value                                                   #
            end                                                                             #
        end                                                                                 #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    for b_index in eachindex(β)                                                             #
        for init_index in eachindex(x_init)                                                 #
            x_before = x_init[init_index]                                                   #
            for time = 1:(size+transient)                                                   #
                x_after = x_before * β[b_index]                                             #
                while (x_after > 1.0)                                                       #
                    x_after = x_after - 1.0                                                 #
                end                                                                         #
                # ------------------------------ TRANSIENT ENDS --------------------------- #
                if (time > transient)                                                       #
                    result[time-transient, 1, init_index, b_index] = x_after                #
                end                                                                         #
                # ------------------------------------------------------------------------- #
                x_before = x_after                                                          #
            end                                                                             #
        end                                                                                 #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    return result, x_init                                                                   #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
#               - Main function                                                             #
function main()                                                                             #
    entropy_data, __desc = generate_betax(series_size, max_entropy_samples)                 #
    net_train_data, net_train_x_init = generate_betax(series_size, net_train_samples;       #
        no_use_random=!net_train_ramdom)                                                    #
    net_test_data, __desc = generate_betax(series_size, net_test_samples;                   #
        check_content=net_train_x_init)                                                     #
    # ------------------------------------------------------------------------------------- #
    #           Info about our process progress                                             #
    index_status = [1, 1]                                                                   #
    save_status(index_status, ARGS[1])                                                      #
    # ------------------------------------------------------------------------------------- #
    #           Save the data =D                                                            #
    save_data(entropy_data, "entropy_data", ARGS[1])                                        #
    save_data(net_train_data, "train_data", ARGS[1])                                        #
    save_data(net_test_data, "test_data", ARGS[1])                                          #
    # ------------------------------------------------------------------------------------- #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
if !(@isdefined(INIT_MARK))
    main()
end