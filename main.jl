"""

Next necessary steps:

- Get DWs in parallel, for all initialization points?

- do directed walks only for the most critical contingency when considering contingencies

- make sure all OPs in the dataset are unique. Get discretization interval 

- make sure all directionaries are only specified in init.jl. If you import .init somewhere, make sure to specify it
in it to keep overview. Also, clean_temp_folder location could be in init. Also DW.txt file in init.

- get another big case to work on (maybe 240 buses), besides 39 and 162


Future work: 

compare method against baselines for transient stability! 
Maybe also voltage collapse; point where you can't solve opf anymore


"""

##########################################################
# This is the main file for the dataset creation toolbox.
##########################################################

# activate current path and initialize environment
using Pkg
using Distributed

# Activate the environment and instantiate packages
Pkg.activate(@__DIR__)
Pkg.instantiate() 
Pkg.status()

# Set the current working directory for all workers (optional but recommended)
cd(@__DIR__)

# Include support scripts on the main process (only if they're not required on workers)
include("functions/method.jl")
include("functions/contingency.jl")
include("functions/acpfcorrect.jl")
include("functions/obbt_lu.jl")
include("functions/polytope.jl")
include("functions/support.jl")
include("functions/write_dfs.jl")
include("functions/directed_walk.jl")
include("functions/ssa_module.jl")
include("functions/dynamics.jl")

# Include initialization module only on the main process
include("init.jl")
using .Initialize

clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

# check if properly initialized
check_initialization()

variable_loads = Initialize.variable_loads

# Record start time
start_time = time()

# Construct polytope using seperating hyperplanes and sample from within the polytope
_, _, _, header = instantiate_system_QCRM(Initialize.network_data, variable_loads) # obtain header of variables
pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1])) # obtain bus numbers of variables
op_header = header_full_op(Initialize.network_data, header, pg_numbers, vm_numbers) # obtain header of full operating points i.e. all pg, vm and pd

# start OBBT, seperating hyperplanes and sampling within polytope
result_init, elapsed_time_init, memory_init, garbage_init = @timed begin
    data_tight_tmp, volumes, hyperplane_time = seperating_hyperplanes(Initialize.network_data, Initialize.hyperplanes, Initialize.variable_loads, Initialize.contingencies_n1, Initialize.stopping_iteration, Initialize.stopping_percentage)
    sample_polytope(Initialize.network_data, data_tight_tmp, Initialize.polytope_samples, Initialize.variable_loads, Initialize.contingencies_n1)
end

feasible_ops_polytope, pf_results_feas_polytope, op_info_feas_pol, infeasible_ops_polytope, pf_results_infeas_polytope, op_info_infeas_pol, nb_feasible, nb_infeasible, pvpq_feasible, initial_feasible, correct_feasible, sampling_time = result_init

# write macros 
df_macros_init(volumes, elapsed_time_init, hyperplane_time, sampling_time, memory_init, initial_feasible, pvpq_feasible, correct_feasible, Initialize.directory, Initialize.pol_macros_filename)

# Perform small-signal stability analysis on all samples
if Initialize.sss_analysis == true

    # define stability boundary
    stability_lower_bound = Initialize.stability_boundary - Initialize.stability_margin
    stability_upper_bound = Initialize.stability_boundary + Initialize.stability_margin

    # perform small signal stability analysis for feasible samples
    result_sss_feas, _, _, _ = @timed begin # elapsed_time_sss, memory_sss, garbage_sss
        sss_evaluation(data_tight_tmp, feasible_ops_polytope, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    damp_pol_feas, dist_pol_feas, _ = result_sss_feas

    # perform small signal stability analysis for infeasible samples
    result_sss_infeas, _, _, _ = @timed begin
        sss_evaluation(data_tight_tmp, infeasible_ops_polytope, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    damp_pol_infeas, dist_pol_infeas, _ = result_sss_infeas
end


if Initialize.directed_walks == true 
    # obtain smalles Pmax to determin minimal distance between OPs R
    min_pmax = 1000
    for i in 1:length(data_tight_tmp["gen"])
        global min_pmax
        pmax = data_tight_tmp["gen"]["$i"]["pmax"]
        if pmax < min_pmax
            min_pmax = pmax
        end
    end

    distance = Initialize.distance  # Distance determining step size
    alpha = Initialize.alpha  # Step size

    # Determine minimal distance between OPs
    R_min = 0.1 # minimum(alpha) * min_pmax * 0.5

    # Get feasible and stable ops with spacing R from each other
    cls_op, closest_ops = remove_nearby_arrays(feasible_ops_polytope, damp_pol_feas, R_min)
    num_dws = length(cls_op)

    if Initialize.dw_computation == "parallel"

        # Perform directed walks in parallel using distributed workers
        start_time_dws = time()

        # Add workers for distributed parallelism
        addprocs(3)  # Adjust based on your needs
        println("Total active workers: ", nworkers())

        # Activate the environment and instantiate packages
        @everywhere begin 
            using Pkg
            Pkg.activate(@__DIR__)
            #@everywhere Pkg.instantiate() 

            # Set the current working directory for all workers (optional but recommended)
            cd(@__DIR__)

            # Include support scripts on the main process (only if they're not required on workers)
            include("functions/method.jl")
            include("functions/contingency.jl")
            include("functions/acpfcorrect.jl")
            include("functions/obbt_lu.jl")
            include("functions/polytope.jl")
            include("functions/support.jl")
            include("functions/write_dfs.jl")
            include("functions/directed_walk.jl")
            include("functions/ssa_module.jl")
            include("functions/dynamics.jl")

            # Include initialization module only on the main process
            include("init.jl")
            using .Initialize

            # clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

        end

        @everywhere begin # $ 'copies'/broadcasts the variable from the local process
            data_tight_tmp = $data_tight_tmp
            feasible_ops_polytope = $feasible_ops_polytope
            closest_ops = $closest_ops
            cls_op = $cls_op
            pf_results_feas_polytope = $pf_results_feas_polytope
            num_dws = length(cls_op)
            variable_loads = $variable_loads
            distance = $distance
            alpha = $alpha

        end

        # Perform directed walks in parallel using Distributed.@distributed
        print("I'm trying to do parallel computing.... fingers crossed!")
        results = @sync @distributed (vcat) for i in 1:num_dws ## TRY WITHOUT @SYNC! or with ASYNC
            dw_ops, dw_stability = DW_step_single_op(data_tight_tmp, feasible_ops_polytope, i, closest_ops[i], cls_op[i], variable_loads, pf_results_feas_polytope, distance, alpha, Initialize.dir_dynamics, Initialize.case_name)
            (dw_ops, dw_stability)  # Return a tuple of results, the last line is collected by the @distributed macro. same as 'return'.
        end

        # remove workers and continue code on one core
        rmprocs(workers())

        # Extracting `dw_ops` and `dw_stability` using comprehension
        directed_walk_ops = vcat([result[1] for result in results]...)
        directed_walk_stability = vcat([result[2] for result in results]...)


        end_time_dws = time()
        elapsed_time_dws = end_time_dws - start_time_dws

    elseif Initialize.dw_computation == "series"

        # Perform directed walks
        start_time_dws = time()
        directed_walk_ops, directed_walk_stability = DW_step(data_tight_tmp, feasible_ops_polytope, closest_ops, cls_op, Initialize.variable_loads, pf_results_feas_polytope, distance, alpha, Initialize.dir_dynamics, Initialize.case_name)
        end_time_dws = time()
        elapsed_time_dws = end_time_dws - start_time_dws

    end

    println("Directed walks completed in $(elapsed_time_dws) seconds")

    # Check AC feasibility of operating points after directed walks
    feasible_ops_dws, pf_results_feas_dws, op_info_feas_dws, infeasible_ops_dws, pf_results_infeas_dws, op_info_infeas_dws, nb_feasible_dws, nb_infeasible_dws, initial_feasible_dws, damp_dws_feas, damp_dws_infeas, dist_dws_feas, dist_dws_infeas = dw_ops_feasibility(Initialize.network_data, data_tight_tmp, Initialize.variable_loads, directed_walk_ops, directed_walk_stability)

end


# Create a multivariate normal distribution and sample from this distribution
if Initialize.mvnd_sampling == true

    if Initialize.sss_analysis == true

        if Initialize.directed_walks == true
            # take AC feasible and small-signal stable samples
            feasible = vcat(feasible_ops_polytope, feasible_ops_dws)
            damp_ops = vcat(damp_pol_feas, damp_dws_feas)
        else
            # take AC feasible and small-signal stable samples
            feasible = vcat(feasible_ops_polytope)
            damp_ops = vcat(damp_pol_feas)
        end

        # get indices of ops in stability boundary region which are AC feasible
        boundary_ops_ind = ops_in_stability_boundary(damp_ops, stability_lower_bound, stability_upper_bound)
        boundary_ops = feasible[boundary_ops_ind]

    else

        boundary_ops = feasible_ops_polytope

    end

    if boundary_ops == []
        throw(ErrorException("There are no AC feasible samples which are small-signal stable in your current dataset. Fiting a MVND is not possible."))
    else
    	# determine number of boundary ops for constructing MVND
    	num_boundary_ops = length(boundary_ops)
    	println("number of boundary ops: ", num_boundary_ops)    
    end

    # start sampling from the multivariate normal distribution, constructed from AC feasible OPs in the stability boundary
    result_mvnd, elapsed_time_mvnd, memory_mvnd, garbage_mvnd = @timed begin
        sample_mvnd(boundary_ops, data_tight_tmp, Initialize.mvnd_samples, Initialize.contingencies_n1) # feasible_ops_polytope
    end

    feasible_ops_mvnd, pf_results_feas_mvnd, op_info_feas_mvnd, infeasible_ops_mvnd, pf_results_infeas_mvnd, op_info_infeas_mvnd, nb_feasible_mvnd, nb_infeasible_mvnd, pvpq_feasible_mvnd, initial_feasible_mvnd, correct_feasible_mvnd, mvnd_sampling_time = result_mvnd

    # write macros
    df_macros_mvnd(Initialize.mvnd_sampling, elapsed_time_mvnd, mvnd_sampling_time, memory_mvnd, initial_feasible_mvnd, pvpq_feasible_mvnd, correct_feasible_mvnd, Initialize.directory, Initialize.mvnd_macros_filename)

    # Perform small-signal stability analysis on all samples
    if Initialize.sss_analysis == true

        # perform small signal stability analysis for both feasible and infeasible samples
        result_sss_mvnd_feas, _, _, _ = @timed begin
            sss_evaluation(data_tight_tmp, feasible_ops_mvnd, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
        end

        damp_mvnd_feas, dist_mvnd_feas, _ = result_sss_mvnd_feas

        # perform small signal stability analysis for both feasible and infeasible samples
        result_sss_mvnd_infeas, _, _, _ = @timed begin
            sss_evaluation(data_tight_tmp, infeasible_ops_mvnd, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
        end

        damp_mvnd_infeas, dist_mvnd_infeas, _ = result_sss_mvnd_infeas

    end
    
end

# Record end time
end_time = time()

# Calculate elapsed time
elapsed_time = end_time - start_time
println("Elapsed time: ", elapsed_time, " seconds")
df_macros_total(elapsed_time, num_boundary_ops, num_dws, elapsed_time_dws, Initialize.directory, Initialize.method_macros)

# write dataframe to CSV
df_DW = construct_df()
output_path_ops = joinpath(Initialize.directory, Initialize.dataset_filename)
CSV.write(output_path_ops, df_DW; delim=';')  

# print summary
summary_result(Initialize.mvnd_sampling, Initialize.directed_walks)

# construct DT training data
df_flow = construct_dt_data()
output_path_dt_data = joinpath(Initialize.directory, Initialize.flows_filename)
CSV.write(output_path_dt_data, df_flow; delim=';')  
