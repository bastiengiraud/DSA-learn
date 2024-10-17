
"""
Importance sampling.

-> I don't know the distribution of the initial input space, assume MVND.

Entropy is a measure of information content. If all the data belongs to stability boundary, entropy is high. You could make a support function computing
the entropy of all the datasets. Given that the security boundary generally falls in the lower probability region of the operating space, a dataspace
containing samples within the boundary region has a high entropy, as entropy is high when the surprise is high (low probability, frequent occuring)

How do you determine the stability boundary? Do you use a bisection method, or do you use LHC and define a set S of samples in the stability boundary?
-> you need a stability boundary for importance sampling. 

In the bisection method, you find samples on both sides of the boundary, and fit a curve between those points. You can use a Gaussian distribution
with a pdf including the distance to the stability boundary. you interpolate between secure and insecure OPs. (used by Florian)
Q: DK1, you need to generate a MVND to sample from, right? in paper, single normal distribution to sample distances? How to obtain OP?
Q: DK2, g(x) is defined by itself. How to obtain g(x)? If you change I, 

you know your boundary 2.75 < damp < 3.25, so you can also sample and choose the points within that boundary. If you do that to construct g(x), you need 
some copulas for exmample to gnerate correlated multivariate distributions to sample data.

Al-Amin uses K-means clustering to segment OPs in clusters. He finds the cluster with the highest entropy, and generates a hypercube around the cluster
with the highest entropy. Then he uses gapsplit to do additional sampling within this cluster.

Alternatively, you fit a MVND over feasible OPs and sample from that. First, do LHC sampling, and classify OPs as secure/insecure and determine damping.
Then, fit a MVND over the OPs falling within the security boundary. (what I'm doing now) (described by Al-Amin as importance sampling)


"""

# activate path and show active packages
using Pkg
root_dir = dirname(dirname(@__DIR__))
Pkg.activate(root_dir)
Pkg.instantiate() 
Pkg.status()

# include support scripts
include(joinpath(root_dir, "functions/method.jl"))
include(joinpath(root_dir, "functions/contingency.jl"))
include(joinpath(root_dir, "functions/ssa_module.jl"))
include(joinpath(root_dir, "functions/dynamics.jl"))
include(joinpath(root_dir, "functions/directed_walk.jl"))
include(joinpath(root_dir, "functions/acpfcorrect.jl"))
include(joinpath(root_dir, "functions/obbt_lu.jl"))
include(joinpath(root_dir, "functions/polytope.jl"))
include(joinpath(root_dir, "functions/support.jl"))
include(joinpath(root_dir, "functions/write_dfs.jl"))

# import initialization module
include(joinpath(root_dir, "init.jl"))
using .Initialize

# check if properly initialized
check_initialization_imp()

variable_loads = Initialize.variable_loads
clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

# Record start time
start_time = time()

# initialize 
ndim, bus_gen, gen_index, min_lim, max_lim = dim_and_limits_variables(Initialize.network_data, Initialize.variable_loads)
pm, N, vars, header = instantiate_system_QCRM(Initialize.network_data, Initialize.variable_loads)
pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))
op_header = header_full_op(Initialize.network_data, header, pg_numbers, vm_numbers)
nb_samples = Initialize.lhc_imp_samples

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

nb_feasible = 0
nb_infeasible = 0 
pvpq_feasible = 0
initial_feasible = 0
correct_feasible = 0
tollerance = 1e-4

data_opf_verif = deepcopy(Initialize.network_data)

for i in 1:nb_samples     

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
        nb_feasible += 1
        initial_feasible += 1
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
        nb_infeasible += 1

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
            #PF_res2 = solve_ac_opf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
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
            nb_feasible += 1
            correct_feasible += 1
            push!(pf_results_feas, PF_res2["solution"]["nw"]["0"])
            push!(load_results_feas, data_opf_verif)
            push!(op_info_feas, op_flag) 
            println("ACPF correct status:", PF_res2["termination_status"] , "\n")
            println("ACPF correct feasibility:", acpfcorrect_feasibility , "\n")
        else 
            # add infeasible initial sample
            nb_infeasible += 1
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

    # perform small signal stability analysis for infeasible samples
    result_sss_infeas, _, _, _ = @timed begin
        sss_evaluation(Initialize.network_data, infeasible_ops_polytope, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    damp_pol_infeas, dist_pol_infeas, _ = result_sss_infeas
end


################# importance sampling ###################

using Statistics 
using Distributions
using LinearAlgebra

nb_imp_samples = Initialize.nb_imp_samples

# define stability boundary
lower_bound = Initialize.stability_boundary - Initialize.stability_margin
upper_bound = Initialize.stability_boundary + Initialize.stability_margin

# get indices of ops in stability boundary region which are AC feasible
boundary_ops_ind = ops_in_stability_boundary(damp_pol_feas, lower_bound, upper_bound)
boundary_ops = feasible_ops[boundary_ops_ind]

# get number of boundary OPs MVND is constructed over
num_boundary_ops = length(boundary_ops)

# distance to stability boundary of current sampled operating point
function boundary_distance(current_op, boundary_op)
    return norm((current_op - boundary_op), 2)
end

function fit_importance_distribution(samples)
    # Calculate mean
    mean_vec = mean(samples, dims=1)
    
    # Calculate covariance matrix
    cov_matrix = cov(samples)

    eps = 1e-6  # Small value to add to the diagonal
    cov_matrix += eps * I

    # Bias the sampling by scaling the covariance matrix
    biased_cov_matrix = 0.25 * cov_matrix

    # Ensure the biased covariance matrix is positive definite
    while !isposdef(biased_cov_matrix)
        eps *= 10
        biased_cov_matrix = bias_factor * (cov_matrix + eps * I)
    end
    
    # Return the new importance distribution (MultivariateNormal)
    return MvNormal(mean_vec[1], biased_cov_matrix)
end


function sample_from_mvn(g_dist::MultivariateNormal, n_samples::Int)
    return rand(g_dist, n_samples)
end

importance_mvn = fit_importance_distribution(boundary_ops) # fit mvnd over boundary samples
importance_samples = sample_from_mvn(importance_mvn, nb_imp_samples)

# check AC feasibility of all samples
pf_results_feas_imp = []
load_results_feas_imp = []
op_info_feas_imp = []

pf_results_infeas_imp = []
load_results_infeas_imp = []
op_info_infeas_imp = []

nb_feasible_imp = 0
nb_infeasible_imp = 0 
pvpq_feasible_imp = 0
initial_feasible_imp = 0
correct_feasible_imp = 0
tollerance = 1e-4

data_opf_verif = deepcopy(Initialize.network_data)

for i in 1:nb_imp_samples   

    for g in eachindex(pg_numbers)
        data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = importance_samples[g,i] 
    end
    for v in eachindex(vm_numbers)
        data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = importance_samples[length(pg_numbers)+v,i] 
    end
    for d in eachindex(pd_numbers)
        data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = importance_samples[length(pg_numbers)+length(vm_numbers)+d,i]  
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
        nb_feasible_imp += 1
        initial_feasible_imp += 1
        push!(pf_results_feas_imp, PF_res0["solution"])
        push!(load_results_feas_imp, data_opf_verif)
        push!(op_info_feas_imp, op_flag) 
        println("initial status:", PF_res0["termination_status"] , "\n")
        println("initial feasibility: ", initial_feasibility)
        #print(initial_feasibility, vm_vio_over, vm_vio_under, sm_vio, "\n")
    else 
        # add infeasible initial sample
        push!(pf_results_infeas_imp, PF_res0["solution"])
        push!(load_results_infeas_imp, data_opf_verif)
        push!(op_info_infeas_imp, op_flag)
        nb_infeasible_imp += 1

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

        x_hat = importance_samples[:,i]
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
            #PF_res2 = solve_ac_opf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
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
            nb_feasible_imp += 1
            correct_feasible_imp += 1
            push!(pf_results_feas_imp, PF_res2["solution"]["nw"]["0"])
            push!(load_results_feas_imp, data_opf_verif)
            push!(op_info_feas_imp, op_flag) 
            println("ACPF correct status:", PF_res2["termination_status"] , "\n")
            println("ACPF correct feasibility:", acpfcorrect_feasibility , "\n")
        else 
            # add infeasible initial sample
            nb_infeasible_imp += 1
            push!(pf_results_infeas_imp, PF_res2["solution"]["nw"]["0"])
            push!(load_results_infeas_imp, data_opf_verif)
            push!(op_info_infeas_imp, op_flag)
            println("ACPF correct status:", PF_res2["termination_status"] , "\n")
            println("ACPF correct feasibility:", acpfcorrect_feasibility , "\n")
        end

        

    end
    

end

print("number of feasible mvn samples: ", nb_feasible_imp, "\n")
    
# get list of feasible operating points x
feasible_ops_imp = []
for i in 1:length(pf_results_feas_imp[:])
    op = get_Full_OP(Initialize.network_data, pf_results_feas_imp[i], load_results_feas_imp[i])
    push!(feasible_ops_imp, op)
end

# get list of infeasible operating points
infeasible_ops_imp = []
for i in 1:length(pf_results_infeas_imp[:])
    op = get_Full_OP(Initialize.network_data, pf_results_infeas_imp[i], load_results_infeas_imp[i])
    push!(infeasible_ops_imp, op)
end


feasible_ops_mvnd = feasible_ops_imp
pf_results_feas_mvnd = pf_results_feas_imp
op_info_feas_mvnd = op_info_feas_imp
infeasible_ops_mvnd = infeasible_ops_imp
pf_results_infeas_mvnd = pf_results_infeas_imp
op_info_infeas_mvnd = op_info_infeas_imp

clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")


# Perform small-signal stability analysis on all samples
if Initialize.sss_analysis == true

    # perform small signal stability analysis for both feasible and infeasible samples
    result_sss_mvnd_feas, _, _, _ = @timed begin
        sss_evaluation(Initialize.network_data, feasible_ops_mvnd, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    damp_mvnd_feas, dist_mvnd_feas, _ = result_sss_mvnd_feas

    # perform small signal stability analysis for both feasible and infeasible samples
    result_sss_mvnd_infeas, _, _, _ = @timed begin
        sss_evaluation(Initialize.network_data, infeasible_ops_mvnd, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    damp_mvnd_infeas, dist_mvnd_infeas, _ = result_sss_mvnd_infeas

end

# # check how many samples fall in stability boundary after importance sampling
# boundary_ops_ind_imp = ops_in_stability_boundary(total_damp_imp, lower_bound, upper_bound)
# boundary_ops_imp = imp_OPs[boundary_ops_ind_imp]

# Record end time
end_time = time()

# Calculate elapsed time
elapsed_time = end_time - start_time
println("Elapsed time: ", elapsed_time, " seconds")
df_macros_total(elapsed_time, num_boundary_ops, 0, 0, Initialize.directory, Initialize.imp_macros)


# write dataframe to CSV
df_DW = construct_df()
output_path_ops = joinpath(Initialize.directory, Initialize.imp_dataset_filename)
CSV.write(output_path_ops, df_DW; delim=';')  

# print summary
summary_result_imp()

# construct DT training data
df_flow = construct_dt_data()
output_path_dt_data = joinpath(Initialize.directory, Initialize.imp_flows_filename)
CSV.write(output_path_dt_data, df_flow; delim=';')  

# clear temp folder
clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

























# for i in 1:10 # n_iters
#     # get current sample index
#     current_op = boundary_ops_ind[i]

#     init_distribution_damp = total_damp[current_op]
#     importance_distribution_damp = 

#     # get pdfs 
#     lhc_pdf = pdf(lhc_mvn)
#     importance_pdf = pdf(importance_distribution)

#     # get samples
#     lhc_samples = sample_from_mvn(lhc_mvn, 1000)
#     importance_samples = sample_from_mvn(importance_mvn, 1000)

#     # get stability indices
#     total_damp, total_dist, total_eigen = result_sss = 
#         sss_evaluation(Initialize.network_data, global_OPs, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)

#     # get importance
#     for i in 1:1000 
#     importance = nothing
#     if boundary_distance(importance_distribution_damp) <= boundary_distance(init_distribution_damp)
#         importance = 1
#     elseif boundary_distance(importance_distribution_damp) > boundary_distance(init_distribution_damp)
#         importance = 0
#     end


# end

