
"""
LHC based sampling, or OPF based sampling. Samples are drawn from input space using LHC sampling. 
An AC PF is solved, if the sample is infeasible, the full AC OPF is solved to find the closest feasible point.

"""

# specify current path
cd("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code")

# activate path and show active packages
using Pkg
Pkg.activate(".")
Pkg.status()

# include support scripts
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/method.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/contingency.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/ssa_module.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/dynamics.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/directed_walk.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/acpfcorrect.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/obbt_lu.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/polytope.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/support.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/write_dfs.jl")

# import initialization module
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/init.jl")
using .Initialize

check_initialization_lhc()
clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

variable_loads = Initialize.variable_loads

# Record start time
start_time = time()

# initialize 
ndim, bus_gen, gen_index, min_lim, max_lim = dim_and_limits_variables(Initialize.network_data, Initialize.variable_loads)
pm, N, vars, header = instantiate_system_QCRM(Initialize.network_data, Initialize.variable_loads)
pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))
op_header = header_full_op(Initialize.network_data, header, pg_numbers, vm_numbers)
nb_samples = Initialize.lhc_samples

# sampling function
function gen_samples_vectors_naive(nb_samples, n_dimensions, level_min, level_max)
    scaling_list = [(level_min[i], level_max[i]) for i in eachindex(level_min)]
    plan, _ = LHCoptim(nb_samples, n_dimensions, 2)
    scaled_plan = scaleLHC(plan, scaling_list)
    return scaled_plan
end

# obtain samples
sample_ops = gen_samples_vectors_naive(nb_samples, ndim, min_lim, max_lim)

# check AC feasibility of all samples
pf_results_feas = []
load_results_feas = []
op_info_feas = []

pf_results_infeas = []
load_results_infeas = []
op_info_infeas = []

global nb_feasible = 0
global nb_infeasible = 0 
global pvpq_feasible = 0
global initial_feasible = 0
global correct_feasible = 0
tollerance = 1e-4

for i in 1:nb_samples    
    data_opf_verif = deepcopy(Initialize.network_data)

    for g in eachindex(pg_numbers)
        data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = sample_ops[i,g] 
    end
    for v in eachindex(vm_numbers)
        data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = sample_ops[i,length(pg_numbers)+v] 
    end
    for d in eachindex(pd_numbers)
        data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = sample_ops[i,length(pg_numbers)+length(vm_numbers)+d]  
        pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
        pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
        sd = pd / pf
        qd = sqrt(sd^2 - pd^2)
        data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
    end

    # dictionary placeholder with OP flags. 1 is feasible, 0 is infeasible
    op_flag = Dict(
        "N0" => 1, # flag for base-case feasibility
        "N1" => 1, # flag for N1 feasibility
        "N1_over_volt" => 0.0, # amount over voltage violation
        "N1_under_volt" => 0.0, # amount under voltage violation
        "N1_flow" => 0.0 # amount flow violation
    )

    # add contingencies
    data_opf_verif["area_gens"] = Dict()
    contingencies = []
    contingencies_gen = []

    # Iterate over the branches in network_data
    for (i, branch) in data_opf_verif["branch"]
        i = parse(Int, i)
        if (i in Initialize.contingencies_n1)
            push!(contingencies, (idx=i, label="LINE-$(i)", type="branch"))
        end
    end

    # add contingencies to tightened network data
    data_opf_verif["branch_contingencies"] = contingencies
    data_opf_verif["gen_contingencies"] = contingencies_gen

    # construct SCOPF in form of multi network formulation
    multinetwork = build_c1_scopf_multinetwork_modif(data_opf_verif)
    # if length(multinetwork["nw"]) == 1
    #     op_flag["N1"] = 0
    # end

    if multinetwork["per_unit"] == true
        for (n, network) in multinetwork["nw"]
            multinetwork["nw"]["$n"]["per_unit"] = true
        end
    end

    # check initial feasibility of base case and contingency cases
    PF_res0 = nothing
    initial_feasibility = nothing
    vm_vio_over = 0.0
    vm_vio_under = 0.0
    sm_vio = 0.0

    for i in 0:(length(multinetwork["nw"])-1)
        # https://github.com/lanl-ansi/PowerModels.jl/blob/30e3d392fd95f1c5a81abd248ac3acc9462211b8/src/prob/pf.jl#L286-L292
        PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)) # solves non linear pf (unbounded), constraints not enforced
        update_data!(multinetwork["nw"]["$i"], PF_res0["solution"]) # update data with PF results
        flows0 = calc_branch_flow_ac(multinetwork["nw"]["$i"]) # compute branch flows
        update_data!(multinetwork["nw"]["$i"], flows0) # add branch flows
        update_data!(PF_res0["solution"], flows0) # add branch flows to solution
        initial_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(Initialize.network_data, PF_res0["solution"], tollerance)

        if i == 0
            # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
            if initial_feasibility != true    
                op_flag["N0"] = 0
            end
        end

        if i != 0
            # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
            if initial_feasibility != true
                op_flag["N1"] = 0
                op_flag["N1_over_volt"] += vm_vio_over
                op_flag["N1_under_volt"] += vm_vio_under
                op_flag["N1_flow"] += sm_vio
            end
        end

    end        

    # check feasibility
    if op_flag["N0"] == 1 && op_flag["N1"] == 1
        global nb_feasible += 1
        global initial_feasible += 1
        push!(pf_results_feas, PF_res0["solution"])
        push!(load_results_feas, data_opf_verif)
        push!(op_info_feas, op_flag) 
        println("initial status:", PF_res0["termination_status"] , "\n")
        println("initial feasibility: ", initial_feasibility)
        #print(initial_feasibility, vm_vio_over, vm_vio_under, sm_vio, "\n")
    else 
        # add infeasible initial sample
        push!(pf_results_infeas, PF_res0["solution"])
        push!(load_results_infeas, data_opf_verif)
        push!(op_info_infeas, op_flag)
        global nb_infeasible += 1

        # op placeholder
        op_flag = Dict(
            "N0" => 1, # flag for base-case feasibility
            "N1" => 1, # flag for N1 feasibility
            "N1_over_volt" => 0.0, # amount over voltage violation
            "N1_under_volt" => 0.0, # amount under voltage violation
            "N1_flow" => 0.0 # amount flow violation
        )

        # if not feasible, map to closest point on feasible region
        vars_new = []

        pm = instantiate_model(multinetwork, ACPPowerModel, build_c1_scopf_load)

        # get generator variables and add to vars
        slack_gen_idx = get_slack_idx_mn(pm)
        for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
            if i != slack_gen_idx && Initialize.network_data["gen"]["$i"]["pmax"] > 0.0
                push!(vars_new, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
            end
        end

        # get voltage magnitude variables and add to vars
        gen_indexes = unique(map(x -> x["gen_bus"], values(pm.data["nw"]["0"]["gen"])))
        for g in gen_indexes
            push!(vars_new, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
        end

        load_vars = collect(variable_loads)
        for d in load_vars
            push!(vars_new, JuMP.variable_by_name(pm.model, string("0_pd[",d,"]")))
        end

        x_hat = sample_ops[i,:]
        N = length(vars_new)

        # additional seperating hyperplanes constraints
        @variable(pm.model, r); # r is the radius of the hypersphere from the sampled point to the relaxation
        @variable(pm.model, aux_new[1:N]) # add auxiliary variable (x_opt)
        @constraint(pm.model, aux_new .== vars_new) # constrain auxiliary variables
        @constraint(pm.model, con_sphere, sqrt(sum((aux_new[i]-x_hat[i])^2 for i in 1:N)) <= r) # add constraint on radius

        @objective(pm.model, Min, r); # the objective is to minimize r

        # solve mapping to feasible region
        PF_res2 = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
        acpfcorrect_feasibility = nothing
        vm_vio_over = 0.0
        vm_vio_under = 0.0
        sm_vio = 0.0

        for i in 0:(length(multinetwork["nw"])-1)
            acpfcorrect_feasibility, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(Initialize.network_data, PF_res2["solution"]["nw"]["$i"], tollerance)

            if i == 0
                # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
                if acpfcorrect_feasibility != true
                    op_flag["N0"] = 0
                end
            end

            if i != 0
                # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
                if acpfcorrect_feasibility != true
                    op_flag["N1"] = 0
                    op_flag["N1_over_volt"] += vm_vio_over
                    op_flag["N1_under_volt"] += vm_vio_under
                    op_flag["N1_flow"] += sm_vio
                end
            end

        end        

        # check feasibility
        if op_flag["N0"] == 1 && op_flag["N1"] == 1
            global nb_feasible += 1
            global correct_feasible += 1
            push!(pf_results_feas, PF_res2["solution"]["nw"]["0"])
            push!(load_results_feas, data_opf_verif)
            push!(op_info_feas, op_flag) 
            println("ACPF correct status:", PF_res2["termination_status"] , "\n")
            println("ACPF correct feasibility:", acpfcorrect_feasibility , "\n")
        else 
            # add infeasible initial sample
            global nb_infeasible += 1
            push!(pf_results_infeas, PF_res2["solution"]["nw"]["0"])
            push!(load_results_infeas, data_opf_verif)
            push!(op_info_infeas, op_flag)
            println("ACPF correct status:", PF_res2["termination_status"] , "\n")
            println("ACPF correct feasibility:", acpfcorrect_feasibility , "\n")
        end

    end

end

print("number of feasible samples: ", nb_feasible, "\n")
    
# get list of feasible operating points x
feasible_ops = []
for i in 1:length(pf_results_feas[:])
    op = get_Full_OP(Initialize.network_data, pf_results_feas[i], load_results_feas[i])
    push!(feasible_ops, op)
end

# get list of infeasible operating points
infeasible_ops = []
for i in 1:length(pf_results_infeas[:])
    op = get_Full_OP(Initialize.network_data, pf_results_infeas[i], load_results_infeas[i])
    push!(infeasible_ops, op)
end


# rename to reuse existing functions
feasible_ops_polytope = feasible_ops
pf_results_feas_polytope = pf_results_feas
op_info_feas_pol = op_info_feas 
infeasible_ops_polytope = infeasible_ops
pf_results_infeas_polytope = pf_results_infeas
op_info_infeas_pol = op_info_infeas

clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")


# Perform small-signal stability analysis on all samples
if Initialize.sss_analysis == true

    # define stability boundary
    stability_lower_bound = Initialize.stability_boundary - Initialize.stability_margin
    stability_upper_bound = Initialize.stability_boundary + Initialize.stability_margin

    # perform small signal stability analysis for feasible samples
    result_sss_feas, _, _, _ = @timed begin # elapsed_time_sss, memory_sss, garbage_sss
        sss_evaluation(Initialize.network_data, feasible_ops_polytope, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    damp_pol_feas, dist_pol_feas, _ = result_sss_feas

    clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

    # perform small signal stability analysis for infeasible samples
    result_sss_infeas, _, _, _ = @timed begin
        sss_evaluation(Initialize.network_data, infeasible_ops_polytope, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    damp_pol_infeas, dist_pol_infeas, _ = result_sss_infeas
end


# Record end time
end_time = time()

# Calculate elapsed time
elapsed_time = end_time - start_time
println("Elapsed time: ", elapsed_time, " seconds")
df_macros_total(elapsed_time, 0, 0, 0, Initialize.directory, Initialize.lhc_macros)

# write dataframe to CSV
df_DW = construct_df()
output_path_ops = joinpath(Initialize.directory, Initialize.lhc_dataset_filename)
CSV.write(output_path_ops, df_DW; delim=';')  

# print summary
summary_result_lhc()

# construct DT training data
df_flow = construct_dt_data()
output_path_dt_data = joinpath(Initialize.directory, Initialize.lhc_flows_filename)
CSV.write(output_path_dt_data, df_flow; delim=';')  

# clear temp folder
clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")
