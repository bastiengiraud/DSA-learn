"""
issues that need fixing

- powersimulationsdynamics.jl does not work on avr IEEET1. Either construct another AVR, or make initialize_avr! function for IEEET1. Now, random AVRTypeI is added which is
similar to IEEET1 AVR (alledgedly).
- powersystems.jl does not have a model for the BPA_GG governor. Now, no governor added.


Next necessary steps:
- think of what to do with dynamic parameters. gaat het gewoon om het laten zien van de toolbox, en verzinnen we wat parameters? Voegen we alsnog de ideal trafo's toe?
willen we het met high penetration of RES doen? Of willen we ook echt test cases die geschikt zijn voor dynamic studies?

- Do we need the N-1 check if we already check for N-1 feasibility during sampling? Modify code so quantative N-1 criteria added during sampling, not afterwards. otherwise double scopf solving.
SCOPF is necessary though for the OPs after the directed walks.
In sample_polytope, de laatste feasibility check is in principe je N-1 check. het oplossen van PFs is hetzelfde als het erna oplossen van een scopf in de N-1 stap. Haal violations
uit sample_polytope gedeelte. Pas mss N-1 aan naar oplossen van PFs in plaats van oplossen scopf.

- when you reach finish point of DW, sample points around it? Check paper -> check points around it, and take smallest step size and steepest ascent.

- Do we compute the damping of the most critical contingency? or only of base case?

- think about how to frame the load uncertainty/forecasting part of the toolbox. cause now its still effectively for one load profile with some load uncertainty. 
if load uncertainty is in the titel, it should have a prominent place in the paper.

- make module ToolBox work


"""

##########################################################
# This is the main file for the dataset creation toolbox.
##########################################################

# specify current path

cd("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code")

# activate path and show active packages
using Pkg
Pkg.activate(".")
Pkg.status()

# include support scripts
include("functions/method.jl")
include("functions/contingency.jl")
include("functions/ssa_module.jl")
include("functions/dynamics.jl")
include("functions/directed_walk.jl")
include("functions/acpfcorrect.jl")
include("functions/obbt_lu.jl")
include("functions/polytope.jl")
include("functions/support.jl")
include("functions/write_dfs.jl")

# import initialization module
include("init.jl")
using .Initialize


# Construct polytope using seperating hyperplanes and sample from within the polytope
_, _, _, header = instantiate_system_QCRM(Initialize.network_data, Initialize.variable_loads) # obtain header of variables
pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1])) # obtain bus numbers of variables
op_header = header_full_op(Initialize.network_data, header, pg_numbers, vm_numbers) # obtain header of full operating points i.e. all pg, vm and pd

# start OBBT, seperating hyperplanes and sampling within polytope
result_init, elapsed_time_init, memory_init, garbage_init = @timed begin
    data_tight_tmp, volumes, hyperplane_time = seperating_hyperplanes(Initialize.network_data, Initialize.hyperplanes, Initialize.variable_loads, Initialize.contingencies_n1, Initialize.stopping_iteration, Initialize.stopping_percentage)
    sample_polytope(Initialize.network_data, data_tight_tmp, Initialize.polytope_samples, Initialize.variable_loads, Initialize.contingencies_n1)
end

feasible_ops_polytope, pf_results_polytope, violation_dict, infeasible_ops_polytope, index_feas_op_polytope, index_infeas_op_polytope, nb_feasible, nb_infeasible, pvpq_feasible, initial_feasible, correct_feasible, sampling_time = result_init

# write macros 
df_macros_init(volumes, elapsed_time_init, hyperplane_time, sampling_time, memory_init, initial_feasible, pvpq_feasible, correct_feasible, Initialize.directory, Initialize.init_macros_filename)

# current global OPs
global_OPs = vcat(feasible_ops_polytope, infeasible_ops_polytope)
index_res = vcat(index_feas_op_polytope, index_infeas_op_polytope)
pf_results_total = vcat(pf_results_polytope)
nb_polytope_ops = length(index_res)

index_feasible = index_feas_op_polytope
index_infeasible = index_infeas_op_polytope

feasible_ops = feasible_ops_polytope
infeasible_ops = infeasible_ops_polytope


# Create a multivariate normal distribution and sample from this distribution
if Initialize.mvnd_sampling == true
    # start sampling from the multivariate normal distribution
    result_mvnd, elapsed_time_mvnd, memory_mvnd, garbage_mvnd = @timed begin
        sample_mvnd(feasible_ops_polytope, Initialize.network_data, data_tight_tmp, Initialize.mvnd_samples, Initialize.contingencies_n1)
    end

    feasible_ops_mvnd, pf_results_mvnd, infeasible_ops_mvnd, index_feas_op_mvnd, index_infeas_op_mvnd, nb_feasible_mvnd, nb_infeasible_mvnd, pvpq_feasible_mvnd, initial_feasible_mvnd, correct_feasible_mvnd, mvnd_sampling_time = result_mvnd

    # write macros
    df_macros_mvnd(Initialize.mvnd_sampling, elapsed_time_mvnd, mvnd_sampling_time, memory_mvnd, initial_feasible_mvnd, pvpq_feasible_mvnd, correct_feasible_mvnd, Initialize.directory, Initialize.mvnd_macros_filename)

    # current global OPs
    global_OPs = vcat(feasible_ops_polytope, infeasible_ops_polytope, feasible_ops_mvnd, infeasible_ops_mvnd)
    index_res = vcat(index_feas_op_polytope, index_infeas_op_polytope, index_feas_op_mvnd.+nb_polytope_ops, index_infeas_op_mvnd.+nb_polytope_ops)
    pf_results_total = vcat(pf_results_polytope, pf_results_mvnd)

    index_feasible = vcat(index_feas_op_polytope, index_feas_op_mvnd.+nb_polytope_ops)
    index_infeasible = vcat(index_infeas_op_polytope, index_infeas_op_mvnd.+nb_polytope_ops)

    feasible_ops = global_OPs[index_feasible]
    infeasible_ops = global_OPs[index_infeasible]
    
end



# Perform a N-1 security check on all feasible samples
if Initialize.contingency_analysis == true 
    # perform N-1 analysis for all AC feasible samples
    result_contingency, elapsed_time_contingency, memory_contingency, garbage_contingency = @timed begin
        N_1_step(Initialize.network_data, feasible_ops, Initialize.contingencies_inf, Initialize.contingencies_n1)
    end

    # obtain sum of violations for all AC feasible samples and all contingency cases
    total_sm_array, total_qg_array, total_over_array, total_under_array, total_pg_array, index_N1_secure = result_contingency  

    sm_index = [mean(arr) for arr in total_sm_array]
    over_index = [mean(arr) for arr in total_over_array]
    under_index = [mean(arr) for arr in total_under_array]

end



# Perform small-signal stability analysis on all samples
if Initialize.sss_analysis == true

    # perform small signal stability analysis for both feasible and infeasible samples
    result_sss, elapsed_time_sss, memory_sss, garbage_sss = @timed begin
        sss_evaluation(data_tight_tmp, global_OPs, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    # obtain the damping ratio, distance to imaginary axis and eigenvalues for all samples
    total_damp, total_dist, total_eigen = result_sss
end



# Perform directed walks on the samples which are closest to the stability boundary
if Initialize.directed_walks == true
    # only do directed walks with OPs closest to stability boundary
    closest_ops = closest_to_zero_indices_no_abs(total_damp, Initialize.eigenvalues_dws)
    cls_op = global_OPs[closest_ops]

    distance = [0.045, 0.02, 0.01] # distance determining step size
    alpha = [0.1, 0.01, 0.1, 0.01] # step size
    directed_walk_ops = DW_step(data_tight_tmp, index_res, closest_ops,  pf_results_total, distance, alpha, Initialize.dir_dynamics, Initialize.case_name)

    # check AC feasibility of operating points after DWs
    feasible_ops_dws, pf_results_dws, infeasible_ops_dws, index_feas_op_dws, index_infeas_op_dws, nb_feasible_dws, nb_infeasible_dws, initial_feasible_dws = dw_ops_feasibility(Initialize.network_data, data_tight_tmp, Initialize.variable_loads, directed_walk_ops, Initialize.contingency_analysis)

    # do N-1 contingency analysis on newly found data points
    if Initialize.contingency_analysis == true
        total_sm_array_dws, total_qg_array_dws, total_over_array_dws, total_under_array_dws, total_pg_array_dws, index_N1_secure_dws = N_1_step(Initialize.network_data, feasible_ops_dws, Initialize.contingencies_inf, Initialize.contingencies_n1)

        sm_index_dws = [mean(arr) for arr in total_sm_array_dws]
        over_index_dws = [mean(arr) for arr in total_over_array_dws]
        under_index_dws = [mean(arr) for arr in total_under_array_dws]
    end

end





# Write all the operating points to a .csv file
# dataframes of global ops
df_DW_f, df_DW_i = dfs_init(feasible_ops, infeasible_ops)

# dataframes including contingency analysis
if Initialize.contingency_analysis == true
    df_DW_f, df_DW_i = dfs_contingency(df_DW_f, df_DW_i, sm_index, over_index, under_index, index_N1_secure)
end

# dataframes including sss analysis
if Initialize.sss_analysis == true
    df_DW_f, df_DW_i = dfs_sss(df_DW_f, df_DW_i, total_damp, total_dist)
end

# dataframes ONLY of dws samples
if Initialize.directed_walks == true
    df_DW_f_dws, df_DW_i_dws = df_dws(df_DW_f, df_DW_i, feasible_ops_dws, infeasible_ops_dws, sm_index_dws, over_index_dws, under_index_dws, index_N1_secure_dws) 
    df_DW_f = vcat(df_DW_f, df_DW_f_dws)
    df_DW_i = vcat(df_DW_i, df_DW_i_dws)
end

# Check if column names feasible&infeasible dataframes match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

# Combine the DataFrames vertically
df_DW = vcat(df_DW_f, df_DW_i)

# write dataframe to CSV
output_path_ops = joinpath(Initialize.directory, Initialize.dataset_filename)
CSV.write(output_path_ops, df_DW)  

# print summary
summary_result(Initialize.mvnd_sampling, Initialize.directed_walks)


