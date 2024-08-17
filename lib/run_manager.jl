#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#       UFPR - 2024
#
#       As the calculation of the maximum entropy for each
#   system takes some time, I create this script to help when
#   I need to stop the process and start it again later.
# ------------------------------------------------------------------
#       ** Required libraries:
#   - DelimitedFiles.jl
# ------------------------------------------------------------------
module RManager
# ------------------------------------------------------------------
#   Libraries:
using DelimitedFiles
# ------------------------------------------------------------------
"""
    Struct in charge of storing the settings of the run.
"""
struct Config
    project::String
    output::String
    tmp::String
end
# ------------------------------------------------------------------
"""
    Struct for a process.
"""
struct Process
    coeffs::Any
    coeffs_index::Vector{Int}
    initial_values::Any
end
# ------------------------------------------------------------------
"""
    Run the main function of a project.
"""
function main(cfg::Config)
    println(" - Running the project '" * cfg.project * "' ... ")
    #---------------------------------------------------------------
    #       Checks for an execution progress status in the temporary 
    #   folder
    print(" - Attempting to recover previous progress ... ")
    rcv = CheckProgress(cfg)
    pcs = Process([], [], [])

    if (!rcv)
        println("\n   - No previous execution was found, creating a new process . . . ")
        pcs = evalfile(cfg.project * "/create_process.jl", [cfg.tmp * "/" * cfg.project])
    end
    println(pcs)
end
# ------------------------------------------------------------------
"""
    Checks for an execution progress status in the temporary folder
"""
function CheckProgress(cfg::Config)
    #---------------------------------------------------------------
    if (!isdir(cfg.tmp))
        return false
    end

    if (!isdir(cfg.tmp * "/" * cfg.project))
        return false
    end

    if (!isfile(cfg.tmp * "/" * cfg.project * "/" * "progress.dat"))
        return false
    end
    #---------------------------------------------------------------
    return true
end
# ------------------------------------------------------------------
"""
"""
function CreateTmpFile(filename::String, content, cfg_tmp::String)
    #---------------------------------------------------------------
    #   Checks if the temporary folder exists.
    if (!isdir(cfg_tmp))
        #   Create the folder . . .
        mkdir(cfg_tmp)
    end
    #---------------------------------------------------------------
    println(content |> typeof)
    println(content |> size)
    open(cfg_tmp * "/" * filename * ".dat", "w") do file
        #   O DelimitedFiles não funciona, vou ter que fazer uma estrutura de arquivo T.T
        writedlm(file, content[:, 1, :, :])
    end
end
# ------------------------------------------------------------------
"""
    Loads the configuration file from 'config.jl'
"""
function load_settings(; path="")
    print(" - Loading settings... ")
    if (isfile(path * "config.jl"))
        include(path * "config.jl")
        println("DONE!")
        return Config(path * project, path * output, path * tmp)
    else
        println("FAIL: Settings not found, repairing...")
        s_cfg = "#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#       UFPR - 2024
#
#       Settings file =D
# ------------------------------------------------------------------
#       Path to the project folder. Note that this folder must have
#   a corresponding main.jl file.
project = \"beta-x/\"
# ------------------------------------------------------------------
#       Path to the folder where the task helper will save the 
#   output files.
output = \"output/\"
# ------------------------------------------------------------------
#       Path to the temporary folder where the task helper will 
#   store the execution progress and associated files.
tmp = \"tmp/\"
# ------------------------------------------------------------------"
        open(path * "config.jl", "w") do file
            write(file, s_cfg)
        end

        return load_settings(; path=path)
    end
end
# ------------------------------------------------------------------
end