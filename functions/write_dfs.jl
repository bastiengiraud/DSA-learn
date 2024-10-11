
using DataFrames
using CSV


function df_macros_total(elapsed_time, num_boundary_ops, num_dws, time_dws, directory, filename)
    #---------- dataframe with macros ----------------
    # make dataframe for run specifications
    df_init = DataFrame(Total_time = elapsed_time)
    df_init.num_mvnd_boundary_ops = num_boundary_ops
    df_init.number_of_dws = num_dws
    df_init.time_of_dws = time_dws

    output_path_macros = joinpath(directory, filename)
    CSV.write(output_path_macros, df_init; delim=';')  

end


function df_macros_init(volumes, elapsed, hp_time, sample_time, memory, init_feasible, pvpq_feasible, correct_feasible, directory, filename)
    #---------- dataframe with macros ----------------
    # make dataframe for run specifications
    df_init = DataFrame(Volumes = volumes)
    df_init.elapsed = fill(elapsed_time_init, nrow(df_init))
    df_init.HP = fill(hyperplane_time, nrow(df_init))
    df_init.sampling = fill(sampling_time, nrow(df_init))
    df_init.memory = fill(memory_init, nrow(df_init))
    df_init.init_feas = fill(initial_feasible, nrow(df_init))
    df_init.pvpq_feas = fill(pvpq_feasible, nrow(df_init))
    df_init.correct_feas = fill(correct_feasible, nrow(df_init))

    output_path_macros = joinpath(directory, filename)
    CSV.write(output_path_macros, df_init; delim=';')  

end


function df_macros_mvnd(mvnd_sampling, elapsed_time_mvnd, mvnd_sampling_time, memory_mvnd, initial_feasible_mvnd, pvpq_feasible_mvnd, correct_feasible_mvnd, directory, mvnd_macros_filename)
    if mvnd_sampling == true
        #---------------- dataframe with macros -------------
        # make dataframe for run specifications
        df_mvnd = DataFrame(
        elapsed = [elapsed_time_mvnd],
        sampling = [mvnd_sampling_time],
        memory = [memory_mvnd],
        init_feas = [initial_feasible_mvnd],
        pvpq_feas = [pvpq_feasible_mvnd],
        correct_feas = [correct_feasible_mvnd])

        output_path_macros_mvnd = joinpath(directory, mvnd_macros_filename)
        CSV.write(output_path_macros_mvnd, df_mvnd; delim=';')  

    end
end


function dfs_init(feasible_ops, infeasible_ops, op_info_feasible, op_info_infeasible)
    #------------- dataframe with ops -----------
    # create dataframe of feasible OPs
    df_DW_f = DataFrame(hcat(feasible_ops...)', Symbol.(op_header[1]))
    df_DW_f.N0 = ones(nrow(df_DW_f))

    # create dataframe of infeasible OPs
    if isempty(infeasible_ops)
        df_DW_i = copy(df_DW_f) 
        empty!(df_DW_i)
    else 
        df_DW_i = DataFrame(hcat(infeasible_ops...)', Symbol.(op_header[1]))
        df_DW_i.N0 = zeros(nrow(df_DW_i))
    end

    flow_viol_feas = []
    flow_viol_infeas = []

    over_volt_feas = []
    over_volt_infeas = []

    under_volt_feas = []
    under_volt_infeas = []

    N1_feas = []
    N1_infeas = []

    for i in 1:length(feasible_ops)
        push!(flow_viol_feas, op_info_feasible[i]["N1_flow"])
        push!(over_volt_feas, op_info_feasible[i]["N1_over_volt"])
        push!(under_volt_feas, op_info_feasible[i]["N1_under_volt"])
        push!(N1_feas, op_info_feasible[i]["N1"])
    end

    for i in 1:length(infeasible_ops)
        push!(flow_viol_infeas, op_info_infeasible[i]["N1_flow"])
        push!(over_volt_infeas, op_info_infeasible[i]["N1_over_volt"])
        push!(under_volt_infeas, op_info_infeasible[i]["N1_under_volt"])
        push!(N1_infeas, op_info_infeasible[i]["N1"])
    end

    # Add N_1 column based on feasibility
    df_DW_f.flow_viol = flow_viol_feas # [1:length(df_DW_f.N0)]
    df_DW_i.flow_viol = flow_viol_infeas # -1 * ones(nrow(df_DW_i))

    df_DW_f.over_volt = over_volt_feas # over_index[1:length(df_DW_f.Feas)]
    df_DW_i.over_volt = over_volt_infeas # -1 * ones(nrow(df_DW_i))

    df_DW_f.under_volt = under_volt_feas # under_index[1:length(df_DW_f.Feas)]
    df_DW_i.under_volt = under_volt_infeas # -1 * ones(nrow(df_DW_i))

    df_DW_f.N1 = N1_feas # under_index[1:length(df_DW_f.Feas)]
    df_DW_i.N1 = N1_infeas # -1 * ones(nrow(df_DW_i))

    return df_DW_f, df_DW_i

end



function dfs_sss(df_DW_f, df_DW_i, total_damp, total_dist)

    nb_feasible = length(df_DW_f.N0)
    nb_total = length(df_DW_f.N0) + length(df_DW_i.N0)

    df_DW_f.damping = total_damp[1:nb_feasible]
    df_DW_f.distance = total_dist[1:nb_feasible]

    df_DW_i.damping = total_damp[(nb_feasible+1):nb_total]
    df_DW_i.distance = total_dist[(nb_feasible+1):nb_total]

    return df_DW_f, df_DW_i
end


function df_dws(df_DW_f, df_DW_i, feasible_ops_dws, infeasible_ops_dws, opf_info_feas_dws, opf_info_infeas_dws)
    # create dataframe of feasible OPs from directed walks
    if isempty(feasible_ops_dws) # create empty dataframe if there are no feasible points
        df_DW_f_dws = copy(df_DW_f) 
        empty!(df_DW_f_dws)
    else
        df_DW_f_dws = DataFrame(hcat(feasible_ops_dws...)', Symbol.(names(df_DW_f)))
    end

    # create dataframe of infeasible OPs from dws
    if isempty(infeasible_ops_dws)
        df_DW_i_dws = copy(df_DW_f_dws) 
        empty!(df_DW_i_dws)
    else 
        df_DW_i_dws = DataFrame(hcat(infeasible_ops_dws...)', Symbol.(names(df_DW_i)))
    end

    flow_viol_feas = []
    flow_viol_infeas = []

    over_volt_feas = []
    over_volt_infeas = []

    under_volt_feas = []
    under_volt_infeas = []

    N1_feas = []
    N1_infeas = []

    for i in 1:length(feasible_ops_dws)
        push!(flow_viol_feas, opf_info_feas_dws[i]["N1_flow"])
        push!(over_volt_feas, opf_info_feas_dws[i]["N1_over_volt"])
        push!(under_volt_feas, opf_info_feas_dws[i]["N1_under_volt"])
        push!(N1_feas, opf_info_feas_dws[i]["N1"])
    end

    for i in 1:length(infeasible_ops_dws)
        push!(flow_viol_infeas, opf_info_infeas_dws[i]["N1_flow"])
        push!(over_volt_infeas, opf_info_infeas_dws[i]["N1_over_volt"])
        push!(under_volt_infeas, opf_info_infeas_dws[i]["N1_under_volt"])
        push!(N1_infeas, opf_info_infeas_dws[i]["N1"])
    end

    # Add N_1 column based on feasibility
    df_DW_f_dws.flow_viol = flow_viol_feas # [1:length(df_DW_f.N0)]
    df_DW_i_dws.flow_viol = flow_viol_infeas # -1 * ones(nrow(df_DW_i))

    df_DW_f_dws.over_volt = over_volt_feas # over_index[1:length(df_DW_f.Feas)]
    df_DW_i_dws.over_volt = over_volt_infeas # -1 * ones(nrow(df_DW_i))

    df_DW_f_dws.under_volt = under_volt_feas # under_index[1:length(df_DW_f.Feas)]
    df_DW_i_dws.under_volt = under_volt_infeas # -1 * ones(nrow(df_DW_i))

    df_DW_f_dws.N1 = N1_feas # under_index[1:length(df_DW_f.Feas)]
    df_DW_i_dws.N1 = N1_infeas # -1 * ones(nrow(df_DW_i))

    return df_DW_f_dws, df_DW_i_dws

end






function construct_df()

    if Initialize.sss_analysis != true && Initialize.directed_walks != true && Initialize.mvnd_sampling == true
        # collect all ops
        feasible_ops = vcat(feasible_ops_polytope, feasible_ops_mvnd)
        infeasible_ops = vcat(infeasible_ops_polytope, infeasible_ops_mvnd)
    
        op_info_feasible = vcat(op_info_feas_pol, op_info_feas_mvnd)
        op_info_infeasible = vcat(op_info_infeas_pol, op_info_infeas_mvnd)
    
        # dataframes of global ops
        df_DW_f, df_DW_i = dfs_init(feasible_ops, infeasible_ops, op_info_feasible, op_info_infeasible)


    elseif Initialize.sss_analysis == true && Initialize.directed_walks != true && Initialize.mvnd_sampling != true
        # collect all ops
        feasible_ops = vcat(feasible_ops_polytope)
        infeasible_ops = vcat(infeasible_ops_polytope)
    
        op_info_feasible = vcat(op_info_feas_pol)
        op_info_infeasible = vcat(op_info_infeas_pol)
    
        # dataframes of global ops
        df_DW_f, df_DW_i = dfs_init(feasible_ops, infeasible_ops, op_info_feasible, op_info_infeasible)
    
        # vcat damping and eigenvalues
        total_damp = vcat(damp_pol_feas, damp_pol_infeas)
        total_dist = vcat(dist_pol_feas, dist_pol_infeas)
    
        # make dataframes
        df_DW_f, df_DW_i = dfs_sss(df_DW_f, df_DW_i, total_damp, total_dist)
    
    
    elseif Initialize.sss_analysis == true && Initialize.directed_walks == true && Initialize.mvnd_sampling != true
        # collect all ops
        feasible_ops = vcat(feasible_ops_polytope, feasible_ops_dws)
        infeasible_ops = vcat(infeasible_ops_polytope, infeasible_ops_dws)
    
        op_info_feasible = vcat(op_info_feas_pol, op_info_feas_dws)
        op_info_infeasible = vcat(op_info_infeas_pol, op_info_infeas_dws)
    
        # dataframes of global ops
        df_DW_f, df_DW_i = dfs_init(feasible_ops, infeasible_ops, op_info_feasible, op_info_infeasible)
    
        # vcat damping and eigenvalues
        total_damp = vcat(damp_pol_feas, damp_dws_feas, damp_pol_infeas, damp_dws_infeas)
        total_dist = vcat(dist_pol_feas, dist_dws_feas, dist_pol_infeas, dist_dws_infeas)
    
        # make dataframes
        df_DW_f, df_DW_i = dfs_sss(df_DW_f, df_DW_i, total_damp, total_dist)
    
    
    # dataframes including sss analysis
    elseif Initialize.sss_analysis == true && Initialize.directed_walks == true && Initialize.mvnd_sampling == true
        # collect all ops
        feasible_ops = vcat(feasible_ops_polytope, feasible_ops_dws, feasible_ops_mvnd)
        infeasible_ops = vcat(infeasible_ops_polytope, infeasible_ops_dws, infeasible_ops_mvnd)
    
        op_info_feasible = vcat(op_info_feas_pol, op_info_feas_dws, op_info_feas_mvnd)
        op_info_infeasible = vcat(op_info_infeas_pol, op_info_infeas_dws, op_info_infeas_mvnd)
    
        # dataframes of global ops
        df_DW_f, df_DW_i = dfs_init(feasible_ops, infeasible_ops, op_info_feasible, op_info_infeasible)
    
        # vcat damping and eigenvalues
        total_damp = vcat(damp_pol_feas, damp_dws_feas, damp_mvnd_feas, damp_pol_infeas, damp_dws_infeas, damp_mvnd_infeas)
        total_dist = vcat(dist_pol_feas, dist_dws_feas, dist_mvnd_feas, dist_pol_infeas, dist_dws_infeas, dist_mvnd_infeas)
    
        # make dataframes
        df_DW_f, df_DW_i = dfs_sss(df_DW_f, df_DW_i, total_damp, total_dist)

    elseif Initialize.sss_analysis == true && Initialize.directed_walks != true && Initialize.mvnd_sampling == true
        # collect all ops
        feasible_ops = vcat(feasible_ops_polytope, feasible_ops_mvnd)
        infeasible_ops = vcat(infeasible_ops_polytope, infeasible_ops_mvnd)
    
        op_info_feasible = vcat(op_info_feas_pol, op_info_feas_mvnd)
        op_info_infeasible = vcat(op_info_infeas_pol, op_info_infeas_mvnd)
    
        # dataframes of global ops
        df_DW_f, df_DW_i = dfs_init(feasible_ops, infeasible_ops, op_info_feasible, op_info_infeasible)
    
        # vcat damping and eigenvalues
        total_damp = vcat(damp_pol_feas, damp_mvnd_feas, damp_pol_infeas, damp_mvnd_infeas)
        total_dist = vcat(dist_pol_feas, dist_mvnd_feas, dist_pol_infeas, dist_mvnd_infeas)
    
        # make dataframes
        df_DW_f, df_DW_i = dfs_sss(df_DW_f, df_DW_i, total_damp, total_dist)
    end
    
    
    # Check if column names feasible&infeasible dataframes match
    if names(df_DW_f) != names(df_DW_i)
        throw(ArgumentError("The DataFrames have different column names"))
    end
    
    # Combine the DataFrames vertically
    df_DW = vcat(df_DW_f, df_DW_i)

    return df_DW

end



function flatten_sample(sample::Dict{String, Any})
    row = Dict{String, Any}()
    for (key, subdict) in sample
        for (subkey, value) in subdict
            # Create new keys combining the outer and inner keys
            row[string(key, "_", subkey)] = value
        end
    end
    return row
end



function construct_dt_data()

    if Initialize.sss_analysis != true && Initialize.directed_walks != true && Initialize.mvnd_sampling != true
        feasible_flow_data = vcat(pf_results_feas_polytope)
        infeasible_flow_data = vcat(pf_results_infeas_polytope)

        feasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in feasible_flow_data]
        infeasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in infeasible_flow_data]
    
        df_flow_f = DataFrame(feasible_flow_data_flat)
        df_flow_i = DataFrame(infeasible_flow_data_flat)

        df_flow_f.feasible = ones(nrow(df_flow_f))
        df_flow_i.feasible = zeros(nrow(df_flow_i))

        df_flow_f.stable = zeros(nrow(df_flow_f))
        df_flow_i.stable = zeros(nrow(df_flow_i))

    elseif Initialize.sss_analysis != true && Initialize.directed_walks != true && Initialize.mvnd_sampling == true
        feasible_flow_data = vcat(pf_results_feas_polytope, pf_results_feas_mvnd)
        infeasible_flow_data = vcat(pf_results_infeas_polytope, pf_results_infeas_mvnd)

        feasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in feasible_flow_data]
        infeasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in infeasible_flow_data]
    
        df_flow_f = DataFrame(feasible_flow_data_flat)
        df_flow_i = DataFrame(infeasible_flow_data_flat)

        df_flow_f.feasible = ones(nrow(df_flow_f))
        df_flow_i.feasible = zeros(nrow(df_flow_i))

        df_flow_f.stable = zeros(nrow(df_flow_f))
        df_flow_i.stable = zeros(nrow(df_flow_i))


    elseif Initialize.sss_analysis == true && Initialize.directed_walks != true && Initialize.mvnd_sampling != true
        feasible_flow_data = vcat(pf_results_feas_polytope)
        infeasible_flow_data = vcat(pf_results_infeas_polytope)

        feasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in feasible_flow_data]
        infeasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in infeasible_flow_data]
    
        df_flow_f = DataFrame(feasible_flow_data_flat)
        df_flow_i = DataFrame(infeasible_flow_data_flat)

        df_flow_f.feasible = ones(nrow(df_flow_f))
        df_flow_i.feasible = zeros(nrow(df_flow_i))
    
        # vcat damping 
        feas_damp = vcat(damp_pol_feas)
        infeas_damp = vcat(damp_pol_infeas)

        # set flag to 1 if the sample has enough damping
        sss_feas_mask = feas_damp .> Initialize.stability_boundary 
        sss_feas_mask = Int.(sss_feas_mask)

        sss_infeas_mask = infeas_damp .> Initialize.stability_boundary 
        sss_infeas_mask = Int.(sss_infeas_mask)

        df_flow_f.stable = sss_feas_mask
        df_flow_i.stable = sss_infeas_mask 


    elseif Initialize.sss_analysis == true && Initialize.directed_walks == true && Initialize.mvnd_sampling != true
        feasible_flow_data = vcat(pf_results_feas_polytope, pf_results_feas_dws)
        infeasible_flow_data = vcat(pf_results_infeas_polytope, pf_results_infeas_dws)

        feasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in feasible_flow_data]
        infeasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in infeasible_flow_data]
    
        df_flow_f = DataFrame(feasible_flow_data_flat)
        df_flow_i = DataFrame(infeasible_flow_data_flat)

        df_flow_f.feasible = ones(nrow(df_flow_f))
        df_flow_i.feasible = zeros(nrow(df_flow_i))
    
        # vcat damping 
        feas_damp = vcat(damp_pol_feas, damp_dws_feas)
        infeas_damp = vcat(damp_pol_infeas, damp_dws_infeas)

        # set flag to 1 if the sample has enough damping
        sss_feas_mask = feas_damp .> Initialize.stability_boundary 
        sss_feas_mask = Int.(sss_feas_mask)

        sss_infeas_mask = infeas_damp .> Initialize.stability_boundary 
        sss_infeas_mask = Int.(sss_infeas_mask)

        df_flow_f.stable = sss_feas_mask
        df_flow_i.stable = sss_infeas_mask 
    
    elseif Initialize.sss_analysis == true && Initialize.directed_walks == true && Initialize.mvnd_sampling == true
        feasible_flow_data = vcat(pf_results_feas_polytope, pf_results_feas_dws, pf_results_feas_mvnd)
        infeasible_flow_data = vcat(pf_results_infeas_polytope, pf_results_infeas_dws, pf_results_infeas_mvnd)

        feasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in feasible_flow_data]
        infeasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in infeasible_flow_data]
    
        df_flow_f = DataFrame(feasible_flow_data_flat)
        df_flow_i = DataFrame(infeasible_flow_data_flat)

        df_flow_f.feasible = ones(nrow(df_flow_f))
        df_flow_i.feasible = zeros(nrow(df_flow_i))
    
        # vcat damping 
        feas_damp = vcat(damp_pol_feas, damp_dws_feas, damp_mvnd_feas)
        infeas_damp = vcat(damp_pol_infeas, damp_dws_infeas, damp_mvnd_infeas)

        # set flag to 1 if the sample has enough damping
        sss_feas_mask = feas_damp .> Initialize.stability_boundary 
        sss_feas_mask = Int.(sss_feas_mask)

        sss_infeas_mask = infeas_damp .> Initialize.stability_boundary 
        sss_infeas_mask = Int.(sss_infeas_mask)

        df_flow_f.stable = sss_feas_mask
        df_flow_i.stable = sss_infeas_mask 

    elseif Initialize.sss_analysis == true && Initialize.directed_walks != true && Initialize.mvnd_sampling == true
        feasible_flow_data = vcat(pf_results_feas_polytope, pf_results_feas_mvnd)
        infeasible_flow_data = vcat(pf_results_infeas_polytope, pf_results_infeas_mvnd)

        feasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in feasible_flow_data]
        infeasible_flow_data_flat = [flatten_sample(sample["branch"]) for sample in infeasible_flow_data]
    
        df_flow_f = DataFrame(feasible_flow_data_flat)
        df_flow_i = DataFrame(infeasible_flow_data_flat)

        df_flow_f.feasible = ones(nrow(df_flow_f))
        df_flow_i.feasible = zeros(nrow(df_flow_i))
    
        # vcat damping 
        feas_damp = vcat(damp_pol_feas, damp_mvnd_feas)
        infeas_damp = vcat(damp_pol_infeas, damp_mvnd_infeas)

        # set flag to 1 if the sample has enough damping
        sss_feas_mask = feas_damp .> Initialize.stability_boundary 
        sss_feas_mask = Int.(sss_feas_mask)

        sss_infeas_mask = infeas_damp .> Initialize.stability_boundary 
        sss_infeas_mask = Int.(sss_infeas_mask)

        df_flow_f.stable = sss_feas_mask
        df_flow_i.stable = sss_infeas_mask 
  
    end

    
    # Check if column names feasible&infeasible dataframes match
    if names(df_flow_f) != names(df_flow_i)
        throw(ArgumentError("The DataFrames have different column names"))
    end
    
    # Combine the DataFrames vertically
    df_flow = vcat(df_flow_f, df_flow_i)

    return df_flow

end






function summary_result(mvnd_sampling, directed_walks)

    println("----------------Summary of samples--------------")
    println("------ Polytope sampling ------")
    println("Number of feasible samples: ", nb_feasible)
    println("Number of infeasible samples: ", nb_infeasible)
    println("Of which initially feasible: ", initial_feasible)
    println("Of which after pvpq feasible: ", pvpq_feasible)
    println("Of which after correction feasible: ", correct_feasible)

    if mvnd_sampling == true
        println("------ MVND sampling ------")
        println("Number of feasible samples: ", nb_feasible_mvnd)
        println("Number of infeasible samples: ", nb_infeasible_mvnd)
        println("Of which initially feasible: ", initial_feasible_mvnd)
        println("Of which after pvpq feasible: ", pvpq_feasible_mvnd)
        println("Of which after correction feasible: ", correct_feasible_mvnd)
    end

    if directed_walks == true
        println("------ DWs ------")
        println("Number of feasible samples: ", nb_feasible_dws)
        println("Number of infeasible samples: ", nb_infeasible_dws)
    end

end



function summary_result_lhc()

    println("----------------Summary of samples--------------")
    println("------ LHC sampling ------")
    println("Number of feasible samples: ", nb_feasible)
    println("Number of infeasible samples: ", nb_infeasible)
    println("Of which initially feasible: ", initial_feasible)
    # println("Of which after pvpq feasible: ", pvpq_feasible)
    println("Of which after correction feasible: ", correct_feasible)

end



function summary_result_opt(mvnd_sampling)

    println("----------------Summary of samples--------------")
    println("------ Optimization based sampling ------")
    println("Number of feasible samples: ", nb_feasible)
    println("Number of infeasible samples: ", nb_infeasible)
    println("Of which initially feasible: ", initial_feasible)
    # println("Of which after pvpq feasible: ", pvpq_feasible)
    println("Of which after correction feasible: ", correct_feasible)

    if mvnd_sampling == true
        println("------ MVND sampling ------")
        println("Number of feasible samples: ", nb_feasible_mvnd)
        println("Number of infeasible samples: ", nb_infeasible_mvnd)
        println("Of which initially feasible: ", initial_feasible_mvnd)
        println("Of which after pvpq feasible: ", pvpq_feasible_mvnd)
        println("Of which after correction feasible: ", correct_feasible_mvnd)
    end

end



function summary_result_imp()

    println("----------------Summary of samples--------------")
    println("------ LHC sampling ------")
    println("Number of feasible samples: ", nb_feasible)
    println("Number of infeasible samples: ", nb_infeasible)
    println("Of which initially feasible: ", initial_feasible)
    # println("Of which after pvpq feasible: ", pvpq_feasible)
    println("Of which after correction feasible: ", correct_feasible)


    println("------ Importance sampling ------")
    println("Number of feasible samples: ", nb_feasible_imp)
    println("Number of infeasible samples: ", nb_infeasible_imp)
    println("Of which initially feasible: ", initial_feasible_imp)
    println("Of which after pvpq feasible: ", pvpq_feasible_imp)
    println("Of which after correction feasible: ", correct_feasible_imp)


end














