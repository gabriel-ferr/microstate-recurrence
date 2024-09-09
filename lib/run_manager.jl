#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado
#                   and Prof. Sérgio Roberto Lopes
#       UFPR - 2024
#
#       As the calculation of the maximum entropy for each
#   system takes some time, I create this script to help when
#   I need to stop the process and start it again later.
# ----------------------------------------------------------------------------------------- #
#       ** Required libraries:                                                              #
#   - JLD2.jl                                                                               #
#   - Dates.jl                                                                              #
# ----------------------------------------------------------------------------------------- #
module RManager                                                                             #
# ----------------------------------------------------------------------------------------- #
#   Libraries:                                                                              #
using JLD2                                                                                  #
using Dates                                                                                 #
# ----------------------------------------------------------------------------------------- #
"""
    Struct in charge of storing the settings of the run.
"""
struct Config                                                                               #
    project::String                                                                         #
    output::String                                                                          #
    tmp::String                                                                             #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
    Run the main function of a project.
"""
function main(cfg::Config)                                                                  #
    #-------------------------------------------------------------------------------------- #
    #       Verify if the project structure, temporary folder and                           #
    #   output folder is okay.                                                              #
    println(" - Verification of the structure in progress ... ")                            #
    #                                                                                       #
    if (!isdir(cfg.project))                                                                #
        println("Couldn't find the project folder!")                                        #
        return                                                                              #
    end                                                                                     #
    #                                                                                       #
    if (!isdir(cfg.tmp))                                                                    #
        println("Couldn't find the temporary folder!")                                      #
        return                                                                              #
    end                                                                                     #
    #                                                                                       #
    if (!isdir(cfg.output))                                                                 #
        println("Couldn't find the output folder!")                                         #
        return                                                                              #
    end                                                                                     #
    #                                                                                       #
    if (!isfile(cfg.project * "/create_process.jl"))                                        #
        println("Couldn't find the create process script!")                                 #
        return                                                                              #
    end                                                                                     #
    #                                                                                       #
    if (!isfile(cfg.project * "/init.jl"))                                                  #
        println("Couldn't find the project initialization script!")                         #
        return                                                                              #
    end                                                                                     #
    #-------------------------------------------------------------------------------------- #
    println(" - Running the project '" * cfg.project * "' ... ")                            #
    #-------------------------------------------------------------------------------------- #
    #       Checks for an execution progress status in the temporary                        #
    #   folder                                                                              #
    println(" - Attempting to recover previous progress ... ")                              #
    #                                                                                       #
    if (isfile(cfg.tmp * "/" * cfg.project * ".project"))                                   #
        if (!isdir(cfg.tmp * "/" * cfg.project))                                            #
            println("   ERROR: The project file was found, but its temporary folder wasn't.")
            println("          Please, remove the project file: `" * cfg.tmp * "/" * cfg.project * ".project`")
            return                                                                          #
        end                                                                                 #
        print("   A previous project was found, do you want to use it or create a new one? (Y - use it | N - create a new one) : ")
        __cmd_status_waiting = true                                                         #
        while __cmd_status_waiting                                                          #
            content = readline()                                                            #
            if (content == "N" || content == "n")                                           #
                rm(cfg.tmp * "/" * cfg.project * ".project")                                #
                rm(cfg.tmp * "/" * cfg.project, recursive=true)                             #
                __cmd_status_waiting = false                                                #
            elseif (content == "Y" || content == "y")                                       #
                __cmd_status_waiting = false                                                #
            else                                                                            #
                print("Please, use `Y` or `N`: ")                                           #
            end                                                                             #
        end                                                                                 #
    end                                                                                     #
    #                                                                                       #
    if (!isfile(cfg.tmp * "/" * cfg.project * ".project"))                                  #
        println("   Creating a new process ...")                                            #
        evalfile(cfg.project * "/create_process.jl", [cfg.tmp * "/" * cfg.project])         #
    end                                                                                     #
    #-------------------------------------------------------------------------------------- #
    #       Run our code =D                                                                 #
    print("\n")                                                                             #
    evalfile(cfg.project * "/init.jl", [cfg.tmp * "/" * cfg.project,                        #
        cfg.output * "/" * cfg.project])                                                    #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
export save_status                                                                          #
"""
        Save the status of a process in the temporary folder.
"""
function save_status(status, path)                                                          #
    save_object(path * ".project", status)                                                  #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
export load_status
"""
        Load status of a process from the temporary folder.
"""
function load_status(path)                                                                  #
    return load_object(path * ".project")                                                   #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
export save_data                                                                            #
"""
        Save a data struct to project temporary folder to be recovery when we need.
"""
function save_data(data, tag, path)                                                         #
    if (!isdir(path))                                                                       #
        mkdir(path)                                                                         #
    end                                                                                     #
    save_object(path * "/" * tag * ".jld2", data)                                           #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
export load_data                                                                            #
"""
        Load a data struct saved in the temporary folder with a tag name.
"""
function load_data(tag, path)                                                               #
    return load_object(path * "/" * tag * ".jld2")                                          #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
export prepare_to_export                                                                    #
"""
        Prepares a folder to export the output data files.
"""
function prepare_to_export(output_path)                                                     #
    if (!isdir(output_path))                                                                #
        mkdir(output_path)                                                                  #
    end                                                                                     #
    path_name = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")                                  #
    mkdir(output_path * "/" * path_name)                                                    #
    return output_path * "/" * path_name                                                    #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
"""
    Loads the configuration file from 'config.jl'
"""
function load_settings(; path="")                                                           #
    print(" - Loading settings... ")                                                        #
    if (isfile(path * "config.jl"))                                                         #
        include(path * "config.jl")                                                         #
        println("DONE!")                                                                    #
        return Config(path * project, path * output, path * tmp)                            #
    else                                                                                    #
        println("FAIL: Settings not found, repairing...")                                   #
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
        open(path * "config.jl", "w") do file                                               #
            write(file, s_cfg)                                                              #
        end                                                                                 #
        #                                                                                   #
        return load_settings(; path=path)                                                   #
    end                                                                                     #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
end