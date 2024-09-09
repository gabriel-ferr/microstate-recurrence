#
#       By Gabriel Ferreira
#       Orientation: Prof. Thiago de Lima Prado 
#                   and Prof. Sérgio Roberto Lopes
#       UFPR - 2024
#
#       This is the actual script that does what we want to do with the Beta-X. 
#   I will first calculate the threshold that gives us the maximum entropy in 
#   the system using different recurrence methods, then we will use it to 
#   calculate the microstate probability of a data set and train an MLP using it
#   to check the difference between the recurrences when we apply them in a 
#   situation like this.
# ----------------------------------------------------------------------------------------- #
#       ** Required libraries:                                                              #
#   - Local: RManager.jl                                                                    #
#   - Local: RP.jl                                                                          #
#   Flux.jl                                                                                 #
#   CairoMakie.jl                                                                           #
#   Statistics.jl                                                                           #
#   LinearAlgebra.jl
#   ProgressMeter.jl                                                                        #
# ----------------------------------------------------------------------------------------- #
INIT_MARK = true                                                                            #
# ----------------------------------------------------------------------------------------- #
if !(@isdefined RManager)                                                                   #
    include("../lib/run_manager.jl")                                                        #
    using .RManager                                                                         #
end                                                                                         #
include("../lib/rp.jl")                                                                     #
include("create_process.jl")                                                                #
# ----------------------------------------------------------------------------------------- #
using .RP                                                                                   #
using Flux                                                                                  #
using CairoMakie                                                                            #
using Statistics                                                                            #
using LinearAlgebra                                                                         #
using ProgressMeter                                                                         #
# ----------------------------------------------------------------------------------------- #
#           Computes the maximum entropy thresholds for all the recurrences we are going    #
#   to compare. The process is saved step by step so that it can be restored if necessary.  #
#   (note: we don't have a backup, only the current step xD)                                #
function calculate_threshold(status)                                                        #
    data = load_data("entropy_data", ARGS[1])                                               #
    std_entropy_threshold = zeros(Float64, length(β))                                       #
    crr_entropy_threshold = zeros(Float64, length(β), 2)                                    #
    # ------------------------------------------------------------------------------------- #
    if (status[2] > 1)                                                                      #
        std_entropy_threshold = load_data("std_threshold", ARGS[1])                         #
        crr_entropy_threshold = load_data("crr_threshold", ARGS[1])                         #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    @showprogress for β_index = status[2]:length(β)                                         #
        std_entropy_threshold[β_index] = RP.threshold_for_max_entroy_std_rrc(data[:, :, :, β_index], 3; no_use_random_tiles=true, macro_accuracy=20, micro_accurary=50)
        crr_entropy_threshold[β_index, :] .= RP.threshold_for_max_entroy_crr_rrc(data[:, :, :, β_index], 3; no_use_random_tiles=true, macro_accuracy=20, micro_accurary=50)
        # --------------------------------------------------------------------------------- #
        save_data(std_entropy_threshold, "std_threshold", ARGS[1])                          #
        save_data(crr_entropy_threshold, "crr_threshold", ARGS[1])                          #
        # --------------------------------------------------------------------------------- #
        status[2] = β_index + 1                                                             #
        save_status(status, ARGS[1])                                                        #
    end                                                                                     #
    status = [2, 1, 1]                                                                      #
    save_status(status, ARGS[1])                                                            #
    compute_data(status)                                                                    #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
function compute_data(status)                                                               #
    train_data = load_data("train_data", ARGS[1])                                           #
    test_data = load_data("test_data", ARGS[1])                                             #
    train_sz = size(train_data)                                                             #
    test_sz = size(test_data)                                                               #
    #                                                                                       #
    powvec = RP.power_vector(3)                                                             #
    # ------------------------------------------------------------------------------------- #
    train_microstates = zeros(Float64, 2^(3 * 3), train_sz[3] * length(β), 2)               #
    test_microstates = zeros(Float64, 2^(3 * 3), test_sz[3] * length(β), 2)                 #
    # ------------------------------------------------------------------------------------- #
    if (status[2] > 1)                                                                      #
        train_microstates = load_data("train_probs", ARGS[1])                               #
        test_microstates = load_data("test_probs", ARGS[1])                                 #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    std_entropy_threshold = load_data("std_threshold", ARGS[1])                             #
    crr_entropy_threshold = load_data("crr_threshold", ARGS[1])                             #
    # ------------------------------------------------------------------------------------- #
    @showprogress for β_index = status[2]:length(β)                                         #
        for sample = status[3]:max(train_sz[3], test_sz[3])                                 #
            # ----------------------------------------------------------------------------- #
            if (sample <= train_sz[3])                                                      #
                # ------------------------------------------------------------------------- #
                #       Calcula os RP                                                       #
                std_rp = RecurrenceMatrix(train_data[:, :, sample, β_index], std_entropy_threshold[β_index])
                crr_rp = RecurrenceMatrix(train_data[:, :, sample, β_index], crr_entropy_threshold[β_index, :]; rrc=RP.crr_recurrence)
                # ------------------------------------------------------------------------- #
                #       Calcula os tiles (probabilidade dos microestados)                   #
                std_probs, _count = MicrostatesProbability(std_rp, 3; no_use_samples=true, powvec=powvec)
                crr_probs, _count = MicrostatesProbability(crr_rp, 3; no_use_samples=true, powvec=powvec)
                #                                                                           #
                train_microstates[collect(keys(std_probs)), sample+(β_index-1)*train_sz[3], 1] .= collect(values(std_probs))
                train_microstates[collect(keys(crr_probs)), sample+(β_index-1)*train_sz[3], 2] .= collect(values(crr_probs))
                # ------------------------------------------------------------------------- #
            end                                                                             #
            # ----------------------------------------------------------------------------- #
            if (sample <= test_sz[3])                                                       #
                # ------------------------------------------------------------------------- #
                #       Calcula os RP                                                       #
                std_rp = RecurrenceMatrix(test_data[:, :, sample, β_index], std_entropy_threshold[β_index])
                crr_rp = RecurrenceMatrix(test_data[:, :, sample, β_index], crr_entropy_threshold[β_index, :]; rrc=RP.crr_recurrence)
                # ------------------------------------------------------------------------- #
                #       Calcula os tiles (probabilidade dos microestados)                   #
                std_probs, _count = MicrostatesProbability(std_rp, 3; no_use_samples=true, powvec=powvec)
                crr_probs, _count = MicrostatesProbability(crr_rp, 3; no_use_samples=true, powvec=powvec)
                #                                                                           #
                test_microstates[collect(keys(std_probs)), sample+(β_index-1)*test_sz[3], 1] .= collect(values(std_probs))
                test_microstates[collect(keys(crr_probs)), sample+(β_index-1)*test_sz[3], 2] .= collect(values(crr_probs))
                # ------------------------------------------------------------------------- #
            end                                                                             #
            # ----------------------------------------------------------------------------- #
            save_data(train_microstates, "train_probs", ARGS[1])                            #
            save_data(test_microstates, "test_probs", ARGS[1])                              #
            status[2] = β_index                                                             #
            status[3] = sample                                                              #
            save_status(status, ARGS[1])                                                    #
            # ----------------------------------------------------------------------------- #
        end                                                                                 #
        status[3] = 1                                                                       #
        save_status(status, ARGS[1])                                                        #
    end                                                                                     #
    # ------------------------------------------------------------------------------------- #
    #           Atualiza o status para a próxima etapa                                      #
    status = [3]                                                                            #
    save_status(status, ARGS[1])                                                            #
    # ------------------------------------------------------------------------------------- #
    finish_project()                                                                        #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
function finish_project()                                                                   #
    #   - Carrega os dados...                                                               #
    train_data = load_data("train_probs", ARGS[1])                                          #
    test_data = load_data("test_probs", ARGS[1])                                            #
    #
    #   - Prepara para salvar os treco                                                      #
    path = prepare_to_export(ARGS[2])                                                       #
    #
    #   - Registro das acurracy...
    accuracy = zeros(Float64, 2, 100)
    #   - Épocas de treinamento.
    epochs = 80
    # ------------------------------------------------------------------------------------- #
    @showprogress Threads.@threads for tries = 1:100                                        #
        # - Faz o treinamento da rede neural e obtém a matrix de confusão resultante...     #
        matrix = use_mlp(train_data, test_data, epochs)                                     #
        # - Prepara o gráfico do uso da recorrência padrão                                  #
        graph = Figure(size=(600, 1200))
        Axis(graph[1, 1], xlabel="Trusty β", ylabel="Predicted β", title="Standard", xticks=1:8, yticks=1:8)
        hm = heatmap!(matrix[:, :, 1], colormap=:BuGn)
        for i = 1:8
            for j = 1:8
                if (matrix[i, j, 1] > 0)
                    text!(i, j, text=string(matrix[i, j, 1]), align=(:center, :center), color=:white)
                end
            end
        end
        # - Gráfico do uso da recorrência de corredor                                       #
        Axis(graph[2, 1], xlabel="Trusty β", ylabel="Predicted β", title="Corridor", xticks=1:8, yticks=1:8)
        hc = heatmap!(matrix[:, :, 2], colormap=:BuGn)
        for i = 1:8
            for j = 1:8
                if (matrix[i, j, 2] > 0)
                    text!(i, j, text=string(matrix[i, j, 2]), align=(:center, :center), color=:white)
                end
            end
        end
        # - Salva os gráficos
        Colorbar(graph[1, 2], hm)
        Colorbar(graph[2, 2], hc)
        save(string(path, "/matrix-", tries, ".png"), graph)
        #
        accuracy[1, tries] = tr(matrix[:, :, 1]) / sum(matrix[:, :, 1])
        accuracy[2, tries] = tr(matrix[:, :, 2]) / sum(matrix[:, :, 2])
    end                                                                                     #
    #
    #       Gráfico da accuracy...
    acc_graph = Figure()
    ax = Axis(acc_graph[1, 1], xlabel="Training", ylabel="Accuracy", title=string("Accuracy of a MLP with ", epochs, " epochs"), xticks=[1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    scatterlines!(accuracy[1, :], label="Standard", marker=:cross, strokewidth=1,
        markersize=10, color=:red, strokecolor=:black)
    scatterlines!(accuracy[2, :], label="Corridor", marker=:cross, strokewidth=1,
        markersize=10, color=:green, strokecolor=:black)
    Legend(acc_graph[2, 1], ax, orientation=:horizontal)
    save(string(path, "/accuracy.png"), acc_graph)
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
function use_mlp(train_data, test_data, epochs)                                             #
    # ------------------------------------------------------------------------------------- #
    #           Prepara as labels que vão ser utilizadas ...                                #
    labels = zeros(Float64, size(train_data)[2])                                            #
    trusty = zeros(Float64, size(test_data)[2])                                             #
    for β_index in eachindex(β)                                                             #
        labels[1+(β_index-1)*net_train_samples:β_index*net_train_samples] .= ones(Float64, net_train_samples) * β[β_index]
    end                                                                                     #
    #                                                                                       #
    for β_index in eachindex(β)                                                             #
        trusty[1+(β_index-1)*net_test_samples:β_index*net_test_samples] .= ones(Float64, net_test_samples) * β[β_index]
    end                                                                                     #
    #                                                                                       #
    labels = Flux.onehotbatch(labels, β)                                                    #
    trusty = Flux.onehotbatch(trusty, β)                                                    #
    # ------------------------------------------------------------------------------------- #
    #           Monta os dados com o Flux para os dois testes.                              #
    std_loader = Flux.DataLoader((train_data[:, :, 1], labels), batchsize=50, shuffle=true) #
    crr_loader = Flux.DataLoader((train_data[:, :, 2], labels), batchsize=50, shuffle=true) #
    # ------------------------------------------------------------------------------------- #
    #           Modelos (são 2, um para cada =V)                                            #
    #       OBS: As duas são iguais
    std_model = Chain(
        Dense(512 => 128, relu),
        Dense(128 => 32, selu),
        Dense(32 => 16, selu),
        Dense(16 => length(β)),
        softmax
    )
    std_model = f64(std_model)
    #                                                                                       #
    crr_model = std_model
    # ------------------------------------------------------------------------------------- #                                                              #
    opt_std_state = Flux.setup(Flux.Adam(0.0001), std_model)
    opt_crr_state = Flux.setup(Flux.Adam(0.0001), crr_model)
    # ------------------------------------------------------------------------------------- #
    std_losses = []
    crr_losses = []
    # ------------------------------------------------------------------------------------- #
    for _ = 1:epochs
        for (x, y) in std_loader
            loss, grads = Flux.withgradient(std_model) do m
                y_hat = m(x)
                Flux.logitcrossentropy(y_hat, y)
            end
            Flux.update!(opt_std_state, std_model, grads[1])
            push!(std_losses, loss)
        end
        #                                                                                   #
        for (x, y) in crr_loader
            loss, grads = Flux.withgradient(crr_model) do m
                y_hat = m(x)
                Flux.logitcrossentropy(y_hat, y)
            end
            Flux.update!(opt_crr_state, crr_model, grads[1])
            push!(crr_losses, loss)
        end
    end
    # ------------------------------------------------------------------------------------- #
    std_result = std_model(test_data[:, :, 1])                                              #
    crr_result = crr_model(test_data[:, :, 2])                                              #
    # ------------------------------------------------------------------------------------- #
    sz = size(std_result)
    res_mat = zeros(Int, sz[1], sz[1], 2)
    #                                                                                       #
    for i = 1:sz[2]                                                                         #
        mx_std = findmax(std_result[:, i])                                                  #
        mx_crr = findmax(crr_result[:, i])                                                  #
        mx_lbl = findmax(trusty[:, i])                                                      #
        #                                                                                   #
        res_mat[mx_lbl[2], mx_std[2], 1] += 1                                               #
        res_mat[mx_lbl[2], mx_crr[2], 2] += 1                                               #
    end                                                                                     #
    #                                                                                       #
    #graph_std = heatmap(res_mat[:, :, 1])
    #graph_crr = heatmap(res_mat[:, :, 2])
    #save("std_res.png", graph_std)
    #save("crr_res.png", graph_crr)

    #println(100 * (tr(res_mat[:, :, 1]) / sum(res_mat[:, :, 1])))
    #println(100 * (tr(res_mat[:, :, 2]) / sum(res_mat[:, :, 2])))
    return res_mat
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
#               - Main function                                                             #
function main()                                                                             #
    status = load_status(ARGS[1])                                                           #
    println(string("Status: ", status))                                                     #
    if (status[1] == 1)                                                                     #
        println("   > Starting script to calculate the threshold for maximum entropy =D")   #
        calculate_threshold(status)                                                         #
    elseif (status[1] == 2)                                                                 #
        compute_data(status)                                                                #
    elseif (status[1] == 3)                                                                 #
        finish_project()                                                                    #
    end                                                                                     #
end                                                                                         #
# ----------------------------------------------------------------------------------------- #
main()