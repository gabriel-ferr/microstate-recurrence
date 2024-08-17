#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#       UFPR - 2024
#
#       This is a main script that uses the run manager. You can
#   use it as a base to create your own script if you need to.
# ------------------------------------------------------------------
#       ** Settings
const β_size = 1000
# ------------------------------------------------------------------
#       This function creates a process and defines the coefficients
#   and initial values used in each step.
function create_process()
    #---------------------------------------------------------------
    #   Values to the β that we want to use.
    β = [2.89, 3.49, 3.99, 4.59, 5.09, 6.19, 6.89, 7.59]
    #---------------------------------------------------------------
    #       Generate default values to calculate ε, which gives the 
    #   maximum entropy.
    #       - initial values to x that will be used to compute the 
    #   max entropy.
    x_init = rand(Float64, 5)
    CreateTmpFile("entropy_data", x_init, ARGS[3])
    #---------------------------------------------------------------
    return Process(β, [0], [""])
end