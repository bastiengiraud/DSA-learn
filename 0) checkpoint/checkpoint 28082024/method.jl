include("acpfcorrect.jl")
include("ssa_module.jl")
include("dynamics.jl")
include("obbt_lu.jl")
include("polytope.jl")
include("support.jl")
include("contingency.jl")

using LatinHypercubeSampling
using JuMP
using Ipopt
using PowerModels
using PowerModelsAnnex
using Suppressor
using Polyhedra
using LinearAlgebra
using CSV
using DataFrames
using PowerModelsSecurityConstrained
using PowerSimulationsDynamics
using PowerSystems
using Statistics 
using Distributions
using DelimitedFiles


PowerModels.silence()
const TOLERANCE = 1e-4;

function instantiate_system_QCRM(data_model, variable_loads)
    pm = instantiate_model(data_model, QCRMPowerModel, build_opf_load) # PowerModels.build_opf)
    load_vars = collect(variable_loads)

    # get number of variables i.e. number of dimension to sample in. In this case, Pg and Pd number of variables
    vars = []
    slack_gen_idx = get_slack_idx(pm)

    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
        if i != slack_gen_idx && data_model["gen"]["$i"]["pmax"] > 0.0
            # print(i, "\n") # i gets printed in ascending order, and the number of gens is the number of gen variables 
            push!(vars, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
        end
    end

    gen_indexes = unique(map(x -> x["gen_bus"], values(pm.data["gen"])))
    for g in gen_indexes # there is a variable for voltage at every bus, but we only want the generator ones, so use gen index
        push!(vars, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
    end


    # add the load vars
    for d in load_vars # there is only a load variable for every load, so again same tactic as for gens
        push!(vars, JuMP.variable_by_name(pm.model, string("0_pd[",d,"]")))
    end

    N = length(vars) # number of variables
    header = get_header_names(vars) # names of every variable

    return pm, N, vars, header
end


function dim_and_limits_variables(network_data, variable_loads)
    ndim = 0
    bus_gen = []
    gen_index = []
    bus_load = []
    load_index = []
    for j in eachindex(1:length(network_data["gen"]))
        push!(bus_gen, network_data["gen"]["$j"]["gen_bus"])
        if "$j" != get_slack_idx(network_data) && network_data["gen"]["$j"]["pmax"] > 0.0 && network_data["gen"]["$j"]["gen_status"] == 1 # if not slack, pmax > 0 and status is active
            push!(gen_index, network_data["gen"]["$j"]["index"])
            ndim += 1
        end
    end
    bus_gen = unique(bus_gen)

    # uncomment for voltage sampling
    ndim += length(bus_gen) # double dimension to account for voltages

    # for j in eachindex(1:length(network_data["load"]))
    #     push!(bus_load, network_data["load"]["$j"]["load_bus"])
    #     if network_data["load"]["$j"]["pd"] > 0.0 
    #         push!(load_index, network_data["load"]["$j"]["index"])
    #         ndim += 1
    #     end
    # end

    load_vars = collect(variable_loads)
    ndim += length(load_vars)

    min_lim = []
    max_lim = []
    for i in gen_index
        push!(min_lim, network_data["gen"]["$i"]["pmin"])
        push!(max_lim, network_data["gen"]["$i"]["pmax"])
    end
    # uncomment for voltage sampling
    for i in bus_gen
        push!(min_lim, network_data["bus"]["$i"]["vmin"])
        push!(max_lim, network_data["bus"]["$i"]["vmax"])
    end
    for i in load_vars #load_index
        push!(min_lim, 0.95*network_data["load"]["$i"]["pd"]) 
        push!(max_lim, 1.05*network_data["load"]["$i"]["pd"]) 
    end
    
    return ndim, bus_gen, gen_index, min_lim, max_lim
end 



function seperating_hyperplanes(network_data , hyperplanes, variable_loads, contingencies_n1, stopping_iteration, stopping_percentage)

    "To change load input space, change in:
    - module_methods.jl -> dim_and_limits_variables
    - acpfcorrect.jl -> variable_load_power_real (for instantiating build_opf_load)"

    #-------------------------- part 1 -------------------------
    # number of bound tightening iteratons
    iter = 3
    volumes = []

    global optimal_setpoints = []
    global original_optimal = []

    # OBBT on voltage magnitude and phase angle difference variables
    # QCLSPowerModel is a strengthened version of the "Quadratic convex" relaxation. An extreme-point encoding of trilinar terms is used along with a constraint to link the lambda variables in multiple trilinar terms
    data_tight_tmp, _, _ = solve_obbt_load_opf!(network_data, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0), max_iter=iter, model_constructor=QCLSPowerModel)
    ndim, _, _, min_lim, max_lim = dim_and_limits_variables(data_tight_tmp, variable_loads)

    # do LHC sampling in tightened space to obtain Nb_HP amount of samples, construct hyperplanes for all samples
    sample_ops = gen_samples_vectors(hyperplanes, ndim, min_lim, max_lim) # size (Nb_HP, x_hat)

    # create a new optimization model pm with tightened bounds, where you add additional constraints to construct seperating hyperplanes
    data_tight_tmp["area_gens"] = Dict()
    contingencies = []
    contingencies_gen = []

    # Iterate over the branches in network_data
    for (i, branch) in data_tight_tmp["branch"]
        i = parse(Int, i)
        if (i in contingencies_n1)
            push!(contingencies, (idx=i, label="LINE-$(i)", type="branch"))
        end
    end

    # add contingencies to tightened network data
    data_tight_tmp["branch_contingencies"] = contingencies
    data_tight_tmp["gen_contingencies"] = contingencies_gen

    # construct SCOPF in form of multi network formulation
    multinetwork = build_c1_scopf_multinetwork_modif(data_tight_tmp)

    if multinetwork["per_unit"] == true
        for (n, network) in multinetwork["nw"]
            multinetwork["nw"]["$n"]["per_unit"] = true
        end
    end

    # instantiate model using QC relaxed formulation of SCOPF with load variation
    pm = instantiate_model(multinetwork, QCRMPowerModel, build_c1_scopf_load) # PowerModels.build_opf)

    # initialize input space with the tightened bounds, polytope = P
    create_scaled_pol(ndim, min_lim , max_lim)
    v = comp_vol()

    println("Initial volume : ", v)
    push!(volumes, v)

    # initialize vector of variables
    vars = []
    global x_hat = []

    # get generator variables and add to vars
    slack_gen_idx = get_slack_idx_mn(pm)
    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
        if i != slack_gen_idx && data_tight_tmp["gen"]["$i"]["pmax"] > 0.0
            push!(vars, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
        end
    end

    # uncomment for voltage sampling
    # get voltage magnitude variables and add to vars
    gen_indexes = (unique(map(x -> x["gen_bus"], values(pm.data["nw"]["0"]["gen"]))))
    for g in gen_indexes
        push!(vars, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
    end

    # added
    load_vars = collect(variable_loads)
    for d in load_vars
        push!(vars, JuMP.variable_by_name(pm.model, string("0_pd[",d,"]")))
    end

    N = length(vars)
    header = get_header_names(vars)
    x_hat = sample_ops[1,:] # to initialize x_hat_p. Instead of rebuilding the optimization model, replace x_hat

    # additional seperating hyperplanes constraints
    @variable(pm.model, r); # r is the radius of the hypersphere from the sampled point to the relaxation
    @variable(pm.model, aux[1:N]) # add auxiliary variable (x_opt)
    @constraint(pm.model, aux .== vars) # constrain auxiliary variables
    @NLparameter(pm.model, x_hat_p[i = 1:N] == x_hat[i]) # updating the parameter has performance gain vs rebuilding the model
    @NLconstraint(pm.model, con_sphere, sqrt(sum((aux[i]-x_hat_p[i])^2 for i in 1:N)) <= r) # add constraint on radius

    @objective(pm.model, Min, r); # the objective is to minimize r

    # initialize volume
    previous_volume = 1
    counter = 0

    hyperplane_time = @elapsed begin
        # construct seperating hyperplanes
        for j in 2:length(sample_ops[:,1]) 
            println("Iteration $(j-1)-----------------------------------")
            result = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)) # solve optimization

            # print_summary(result["solution"])
        
            if result["termination_status"] == LOCALLY_SOLVED # check if solved
                x_opt = JuMP.value.(vars)
                r_opt = JuMP.value(r)
                push!(optimal_setpoints, x_opt)

                println("x_hat: ", x_hat) # x_hat is the drawn sample
                println("x_opt: ", x_opt) # x_opt is the solution to the optimization
                println("hypershpere radius i.e. QC infeasible: ", r_opt)
            else
                print("x_opt not found. Termination status: ", result["termination_status"], "\n")
                r_opt = 0
            end

            if !(isapprox(JuMP.value(r), 0; atol=TOLERANCE)) && (r_opt > 0) # if optimization solved and initial point was infeasible
                # get transposed normal vector
                n_normT = transpose(x_hat - x_opt[:] )
                normal_normalized = norm(n_normT) # normalize hyperplane to prevent numerical instabilities
                A = n_normT/normal_normalized
                b = n_normT*x_opt/normal_normalized

                # Update results
                add_to_pol(A, b)
                v = comp_vol()

                if v == -1.0
                    print("Invalid hyperplane. Removing ...")
                    remove_from_pol()
                    v = comp_vol()
                end
                
                println("Volume : ", v)
                push!(volumes, v)
            else
                push!(original_optimal, x_hat)

            end

            global x_hat = sample_ops[j,:] # set x_hat for next iteration
            JuMP.set_value.(x_hat_p, x_hat)

            # compute percentage decrease in volume
            decrease_percentage = (previous_volume - v) / previous_volume

            # Update the previous volume to the new volume
            previous_volume = v

            # Check if the decrease is less than 10%
            if decrease_percentage < stopping_percentage
                counter += 1
                println("Volume decrease of less than $(stopping_percentage*100)%. Current consecutive: ", counter, "\n")
            else
                counter = 0  # Reset counter if the decrease is more than 10%
                println(" ")
            end

            # Break the loop if the counter reaches 10
            if counter >= stopping_iteration
                println("Stopping early as the volume decrease is less than $(stopping_percentage*100)% for $stopping_iteration consecutive iterations.")
                break
            end

        end
    end


    return data_tight_tmp, volumes, hyperplane_time

end

function sample_polytope(network_data, data_tight_tmp, polytope_samples, variable_loads, contingencies_n1)

    nb_samples = polytope_samples
    global in_samples = sample_pol(nb_samples) # draw samples from polytope P
    global nb_feasible = 0
    global nb_infeasible = 0
    global pvpq_feasible = 0
    global initial_feasible = 0
    global correct_feasible = 0
    tollerance = 1e-4

    pm, N, vars, header = instantiate_system_QCRM(data_tight_tmp, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))
    pf_results = []
    load_results = []

    sampling_time = @elapsed begin
        for i in 1:nb_samples # loop over number of samples in polytope
            #println("$i______________________")
            data_opf_verif = deepcopy(data_tight_tmp) # copy tightened bounds network data, replace values with drawn samples
            
            # replace values in dictionary with samples drawn from polytope
            for g in eachindex(pg_numbers)
                data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = in_samples[g,i] 
            end
            for v in eachindex(vm_numbers)
                data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = in_samples[length(pg_numbers)+v,i] 
            end
            for d in eachindex(pd_numbers)
                data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = in_samples[length(pg_numbers)+length(vm_numbers)+d,i]  
                pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
                pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
                sd = pd / pf
                qd = sqrt(sd^2 - pd^2)
                data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
            end

            print("Current polytope sample number: ", i, "\n")

            # construct SCOPF in form of multi network formulation
            multinetwork = build_c1_scopf_multinetwork_modif(data_opf_verif)

            if multinetwork["per_unit"] == true
                for (n, network) in multinetwork["nw"]
                    multinetwork["nw"]["$n"]["per_unit"] = true
                end
            end

            # check initial feasibility of base case and contingency cases
            init_feas = 0
            PF_res0 = nothing
            for i in 0:(length(multinetwork["nw"])-1)
                PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
                initial_feasibility = check_ac_feasibility(data_tight_tmp, PF_res0, tollerance)
                if PF_res0["termination_status"] == LOCALLY_SOLVED && initial_feasibility != true
                    init_feas += 1
                end
            end                

            # check feasibility
            if init_feas == 0 
                global nb_feasible += 1
                global initial_feasible += 1
                push!(pf_results, PF_res0["solution"])
                push!(load_results, data_opf_verif)
                println("initial status:", PF_res0["termination_status"] , "\n")
            else 
                # add infeasible initial sample
                push!(pf_results, PF_res0["solution"])
                push!(load_results, data_opf_verif)
                global nb_infeasible += 1

                # solve pvpq
                pvpq_feas = 0
                PF_res1 = nothing
                for i in 0:(length(multinetwork["nw"])-1)
                    PF_res1 = adjust_PVPQ(multinetwork["nw"]["$i"], 4)
                    PVPQ_feasibility = check_ac_feasibility(data_tight_tmp, PF_res1, tollerance)
                    if PF_res1["termination_status"] == LOCALLY_SOLVED && PVPQ_feasibility != true
                        pvpq_feas += 1
                    end
                end

                if pvpq_feas == 0 
                    # if pvpq feasible add feasible point
                    global nb_feasible += 1
                    global pvpq_feasible += 1
                    push!(pf_results, PF_res1["solution"])
                    push!(load_results, data_opf_verif)
                    println("PV/PQ status:", PF_res1["termination_status"] , "\n")
                else
                    # add infeasible pvpq
                    push!(pf_results, PF_res1["solution"])
                    push!(load_results, data_opf_verif)
                    global nb_infeasible += 1

                    vars_new = []

                    pm = instantiate_model(multinetwork, ACPPowerModel, build_c1_scopf_load)

                    # get generator variables and add to vars
                    slack_gen_idx = get_slack_idx_mn(pm)
                    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
                        if i != slack_gen_idx && data_tight_tmp["gen"]["$i"]["pmax"] > 0.0
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

                    x_hat = in_samples[:,i]
                    N = length(vars_new)

                    # additional seperating hyperplanes constraints
                    @variable(pm.model, r); # r is the radius of the hypersphere from the sampled point to the relaxation
                    @variable(pm.model, aux_new[1:N]) # add auxiliary variable (x_opt)
                    @constraint(pm.model, aux_new .== vars_new) # constrain auxiliary variables
                    @NLconstraint(pm.model, con_sphere, sqrt(sum((aux_new[i]-x_hat[i])^2 for i in 1:N)) <= r) # add constraint on radius

                    @objective(pm.model, Min, r); # the objective is to minimize r

                    PF_res2 = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
                    acpfcorrect_feasibility = check_ac_feasibility_mn(data_tight_tmp, PF_res2, tollerance)

                    if PF_res2["termination_status"] == LOCALLY_SOLVED || acpfcorrect_feasibility == true 
                        global nb_feasible += 1
                        global correct_feasible += 1
                        push!(pf_results, PF_res2["solution"]["nw"]["0"])
                        push!(load_results, data_opf_verif)
                        print("acpfcorrect feasibility: ", acpfcorrect_feasibility, "\n")
                        println("ACPF correct status:", PF_res2["termination_status"] , "\n")
                    else
                        push!(pf_results, PF_res2["solution"]["nw"]["0"])
                        push!(load_results, data_opf_verif)
                        global nb_infeasible += 1
                        print("infeasible acpfcorrect added", "\n")
                        println("PV/PQ status:", PF_res1["termination_status"] )
                        println("ACPF correct status:", PF_res2["termination_status"] )
                    end
                end
            end
    
        end
    end

    global over_array = []
    global under_array = []
    global pg_array = []
    global qg_array = []
    global sm_array = []

    # check what kind of violations each sample has
    for i in 1:length(pf_results[:]) 
        # check for voltage violations 
        vm_vio_over, vm_vio_under = check_vm_violations(data_tight_tmp, pf_results[i], tollerance)

        # check for generator violations
        pg_vio, qg_vio = check_pg_pq_violations(data_tight_tmp, pf_results[i], tollerance)

        # check for line violations
        sm_vio = check_flow_violations(data_tight_tmp, pf_results[i], tollerance)
            
        push!(over_array, vm_vio_over) 
        push!(under_array, vm_vio_under) 

        push!(pg_array, pg_vio)
        push!(qg_array, qg_vio)

        push!(sm_array, sm_vio)
    end

    violation_dict = Dict(
        "over_array" => over_array,
        "under_array" => under_array,
        "pg_array" => pg_array,
        "qg_array" => qg_array,
        "sm_array" => sm_array
    )

    print("number of feasible samples: ", nb_feasible, "\n")
    
    # count zeros and get index
    _, index_pg = find_zeros_and_indexes(pg_array, tollerance)
    _, index_qg = find_zeros_and_indexes(qg_array, tollerance)
    _, index_sm = find_zeros_and_indexes(sm_array, tollerance)
    _, index_ovi = find_zeros_and_indexes(over_array, tollerance)
    _, index_uvi = find_zeros_and_indexes(under_array, tollerance)

    # if all indices are zero i.e. no violation
    index_feas_op_polytope = [x for x in index_qg if (x in index_sm  && x in index_ovi && x in index_uvi && x in index_pg)]

    # if there is a violation
    index_infeas_op_polytope = [x for x in eachindex(pf_results) if !(x in index_feas_op_polytope)]

    # get list of feasible operating points x
    feasible_ops_polytope = []
    for i in index_feas_op_polytope
        op = get_Full_OP(data_tight_tmp, pf_results[i], load_results[i])
        push!(feasible_ops_polytope, op)
    end

    # get list of infeasible operating points
    infeasible_ops_polytope = []
    for i in index_infeas_op_polytope
        op = get_Full_OP(data_tight_tmp, pf_results[i], load_results[i])
        push!(infeasible_ops_polytope, op)
    end

    return feasible_ops_polytope, pf_results, violation_dict, infeasible_ops_polytope, index_feas_op_polytope, index_infeas_op_polytope, nb_feasible, nb_infeasible, pvpq_feasible, initial_feasible, correct_feasible, sampling_time
end



function sample_mvnd(feasible_ops_polytope, network_basic, data_tight_tmp, mvnd_samples, contingencies_n1)
    mean_vec = mean(feasible_ops_polytope, dims=1)

    cov_matrix = cov(feasible_ops_polytope)
    eps = 1e-6  # Small value to add to the diagonal
    cov_matrix += eps * I

    # Bias the sampling by scaling the covariance matrix
    biased_cov_matrix = 0.25 * cov_matrix

    # Ensure the biased covariance matrix is positive definite
    while !isposdef(biased_cov_matrix)
        eps *= 10
        biased_cov_matrix = bias_factor * (cov_matrix + eps * I)
    end

    # Create a Multivariate Normal Distribution with the biased covariance matrix
    mvn = MvNormal(mean_vec[1], cov_matrix)

    # Draw N3 samples from the biased distribution
    biased_samples = rand(mvn, mvnd_samples)
    pm, N, vars, header = instantiate_system_QCRM(data_tight_tmp, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))
    pf_results = []
    load_results = []
    
    global nb_feasible_mvnd = 0
    global nb_infeasible_mvnd = 0
    global initial_feasible_mvnd = 0
    global pvpq_feasible_mvnd = 0
    global correct_feasible_mvnd = 0
    tollerance = 1e-4

    mvnd_sampling_time = @elapsed begin
        for i in eachindex(biased_samples[1,:])
            data_opf_verif = deepcopy(data_tight_tmp)
            
            for g in eachindex(pg_numbers)
                data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = biased_samples[g,i] #optimal_setpoints[i][g] #in_samples[g,i]
            end
            for v in eachindex(vm_numbers)
                data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = biased_samples[length(pg_numbers)+v,i] #optimal_setpoints[i][length(pg_numbers)+v] #in_samples[length(pg_numbers)+v,i]
            end
            for d in eachindex(pd_numbers)
                var_load_index = pd_numbers[d]
                data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = biased_samples[length(pg_numbers)+length(vm_numbers)+var_load_index,i] 
                pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
                pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
                sd = pd / pf
                qd = sqrt(sd^2 - pd^2)
                data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
            end

            print("Current polytope sample number: ", i, "\n")

            # construct SCOPF in form of multi network formulation
            multinetwork = build_c1_scopf_multinetwork_modif(data_opf_verif)

            if multinetwork["per_unit"] == true
                for (n, network) in multinetwork["nw"]
                    multinetwork["nw"]["$n"]["per_unit"] = true
                end
            end

            # check initial feasibility of base case and contingency cases
            init_feas = 0
            PF_res0 = nothing
            for i in 0:(length(multinetwork["nw"])-1)
                PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
                initial_feasibility = check_ac_feasibility(data_tight_tmp, PF_res0, tollerance)
                if PF_res0["termination_status"] == LOCALLY_SOLVED && initial_feasibility != true
                    init_feas += 1
                end
            end                

            # check feasibility
            if init_feas == 0 
                global nb_feasible_mvnd += 1
                global initial_feasible_mvnd += 1
                push!(pf_results, PF_res0["solution"])
                push!(load_results, data_opf_verif)
                println("initial status:", PF_res0["termination_status"] , "\n")
            else 
                # add infeasible initial sample
                push!(pf_results, PF_res0["solution"])
                push!(load_results, data_opf_verif)
                global nb_infeasible_mvnd += 1

                # solve pvpq
                pvpq_feas = 0
                PF_res1 = nothing
                for i in 0:(length(multinetwork["nw"])-1)
                    PF_res1 = adjust_PVPQ(multinetwork["nw"]["$i"], 4)
                    PVPQ_feasibility = check_ac_feasibility(data_tight_tmp, PF_res1, tollerance)
                    if PF_res1["termination_status"] == LOCALLY_SOLVED && PVPQ_feasibility != true
                        pvpq_feas += 1
                    end
                end

                if pvpq_feas == 0 
                    # if pvpq feasible add feasible point
                    global nb_feasible_mvnd += 1
                    global pvpq_feasible_mvnd += 1
                    push!(pf_results, PF_res1["solution"])
                    push!(load_results, data_opf_verif)
                    println("PV/PQ status:", PF_res1["termination_status"] , "\n")
                else
                    # add infeasible pvpq
                    push!(pf_results, PF_res1["solution"])
                    push!(load_results, data_opf_verif)
                    global nb_infeasible_mvnd += 1

                    vars_new = []

                    pm = instantiate_model(multinetwork, ACPPowerModel, build_c1_scopf_load)

                    # get generator variables and add to vars
                    slack_gen_idx = get_slack_idx_mn(pm)
                    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
                        if i != slack_gen_idx && data_tight_tmp["gen"]["$i"]["pmax"] > 0.0
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

                    x_hat = in_samples[:,i]
                    N = length(vars_new)

                    # additional seperating hyperplanes constraints
                    @variable(pm.model, r); # r is the radius of the hypersphere from the sampled point to the relaxation
                    @variable(pm.model, aux_new[1:N]) # add auxiliary variable (x_opt)
                    @constraint(pm.model, aux_new .== vars_new) # constrain auxiliary variables
                    @NLconstraint(pm.model, con_sphere, sqrt(sum((aux_new[i]-x_hat[i])^2 for i in 1:N)) <= r) # add constraint on radius

                    @objective(pm.model, Min, r); # the objective is to minimize r

                    PF_res2 = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
                    acpfcorrect_feasibility = check_ac_feasibility_mn(data_tight_tmp, PF_res2, tollerance)

                    if PF_res2["termination_status"] == LOCALLY_SOLVED || acpfcorrect_feasibility == true 
                        global nb_feasible_mvnd += 1
                        global correct_feasible_mvnd += 1
                        push!(pf_results, PF_res2["solution"]["nw"]["0"])
                        push!(load_results, data_opf_verif)
                        print("acpfcorrect feasibility: ", acpfcorrect_feasibility, "\n")
                        println("ACPF correct status:", PF_res2["termination_status"] , "\n")
                    else
                        push!(pf_results, PF_res2["solution"]["nw"]["0"])
                        push!(load_results, data_opf_verif)
                        global nb_infeasible_mvnd += 1
                        print("infeasible acpfcorrect added", "\n")
                        println("PV/PQ status:", PF_res1["termination_status"] )
                        println("ACPF correct status:", PF_res2["termination_status"] )
                    end
                end
            end
        end


    end

    global over_array = []
    global under_array = []
    global pg_array = []
    global qg_array = []
    global sm_array = []
    global res_arr = []

    # check what kind of violations each sample has
    for i in 1:length(pf_results[:]) 
        # check for voltage violations 
        vm_vio_over, vm_vio_under = check_vm_violations(data_tight_tmp, pf_results[i], tollerance)

        # check for generator violations
        pg_vio, qg_vio = check_pg_pq_violations(data_tight_tmp, pf_results[i], tollerance)

        # check for line violations
        sm_vio = check_flow_violations(data_tight_tmp, pf_results[i], tollerance)
            
        push!(over_array, vm_vio_over) 
        push!(under_array, vm_vio_under) 

        push!(pg_array, pg_vio)
        push!(qg_array, qg_vio)

        push!(sm_array, sm_vio)
    end

    nb_pg, index_pg = find_zeros_and_indexes(pg_array, tollerance)
    nb_qg, index_qg = find_zeros_and_indexes(qg_array, tollerance)
    nb_sm, index_sm = find_zeros_and_indexes(sm_array, tollerance)
    nb_ovi, index_ovi = find_zeros_and_indexes(over_array, tollerance)
    nb_uvi, index_uvi = find_zeros_and_indexes(under_array, tollerance)

    print("number of feasible mvnd samples: ", nb_feasible_mvnd, "\n")

    index_feas_op_mvnd = [x for x in index_qg if (x in index_sm  && x in index_ovi && x in index_uvi && x in index_pg)]
    index_infeas_op_mvnd = [x for x in eachindex(pf_results) if !(x in index_feas_op_mvnd)]

    feasible_ops_mvnd = []
    for i in index_feas_op_mvnd
        op = get_Full_OP(data_tight_tmp, pf_results[i], load_results[i])
        push!(feasible_ops_mvnd, op)
    end

    infeasible_ops_mvnd = []
    for i in index_infeas_op_mvnd
        op = get_Full_OP(data_tight_tmp, pf_results[i], load_results[i])
        push!(infeasible_ops_mvnd, op)
    end

    return feasible_ops_mvnd, pf_results, infeasible_ops_mvnd, index_feas_op_mvnd, index_infeas_op_mvnd, nb_feasible_mvnd, nb_infeasible_mvnd, pvpq_feasible_mvnd, initial_feasible_mvnd, correct_feasible_mvnd, mvnd_sampling_time

end





function N_1_step(network_data, OPs_Feas, cont_out, contingencies_n1)

    global over_array = []
    global under_array = []
    global pg_array = []
    global qg_array = []
    global sm_array = []
    global res_arr = []

    global total_over_array = []
    global total_under_array = []
    global total_pg_array = []
    global total_qg_array = []
    global total_sm_array = []

    global N_1_feasible = []

    threshold = 1e-4
    tollerance = 1e-4

    # get number of variables
    pm, N, vars, header = instantiate_system_QCRM(network_data, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))

    data_opf_verif = deepcopy(network_data)

    for i in eachindex(OPs_Feas) # loop over every operating point
        OP = OPs_Feas[i] 
        data_opf_verif = deepcopy(data_opf_verif)
        
        for g in eachindex(pg_numbers)
            data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = OP[g] # set generator setpoints to operating point
        end
        for v in eachindex(vm_numbers)
            data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = OP[length(pg_numbers)+v] # set voltages to operating point
        end
        for d in eachindex(pd_numbers)
            var_load_index = pd_numbers[d]
            data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = OP[length(pg_numbers)+length(vm_numbers)+var_load_index] 
            pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
            pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
            sd = pd / pf
            qd = sqrt(sd^2 - pd^2)
            data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
        end

        print("Current contingency sample number: ", i, "\n")

        global network_data = deepcopy(data_opf_verif)

        network_data["area_gens"] = Dict()
        network_data["gen_contingencies"] = []
        network_data["branch_contingencies"] = []
        contingencies = []
        contingencies_gen = []

        # Iterate over the branches in network_data
        for (i, branch) in network_data["branch"]
            i = parse(Int, i)

            # Check if contingencies_n1 is empty
            if isempty(contingencies_n1)
                # If contingencies_n1 is empty, include all branches not in cont_out
                if !(i in cont_out)
                    push!(contingencies, (idx=i, label="LINE-$(i)", type="branch"))
                end
            else
                # If contingencies_n1 is not empty, include contingencies in contingencies_n1 except for those in cont_out
                if !(i in cont_out) && (i in contingencies_n1)
                    push!(contingencies, (idx=i, label="LINE-$(i)", type="branch"))
                end
            end
        end

        network_data["branch_contingencies"] = contingencies
        network_data["gen_contingencies"] = contingencies_gen

        global multinetwork = build_c1_scopf_multinetwork_modif(network_data)

        if multinetwork["per_unit"] == true
            for (n, network) in multinetwork["nw"]
                multinetwork["nw"]["$n"]["per_unit"] = true
            end
        end
        
        result = run_c1_scopf_load(multinetwork, ACPPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
        for l in 0:(length(multinetwork["nw"])-1)
            update_data!(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
        end    
        
        global result = run_c1_scopf_load(multinetwork, ACPPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
        

        for l in 0:(length(multinetwork["nw"])-1)
            pg_vio, qg_vio = check_pg_pq_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"], tollerance)
            push!(pg_array, pg_vio)
            push!(qg_array, qg_vio)

            sm_vio = check_flow_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"], tollerance)
            push!(sm_array, sm_vio)

            vm_over_tmp, vm_under_tmp = check_vm_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"], tollerance)
            push!(over_array, vm_over_tmp)
            push!(under_array, vm_under_tmp)
            
        end

        push!(total_over_array,over_array)
        push!(total_under_array,under_array)
        push!(total_pg_array,pg_array)
        push!(total_qg_array, qg_array)
        push!(total_sm_array, sm_array)

        if (all_below_threshold(over_array, threshold) &&
            all_below_threshold(under_array, threshold) &&
            #all_below_threshold(pg_array, threshold) &&
            #all_below_threshold(qg_array, threshold) &&
            all_below_threshold(sm_array, threshold)
            )
            push!(N_1_feasible, i)
        end

        
        # Reset local arrays
        over_array = []
        under_array = []
        pg_array = []
        qg_array = []
        sm_array = []

    end
    return total_sm_array, total_qg_array, total_over_array, total_under_array , total_pg_array, N_1_feasible #, result["solution"]["nw"], multinetwork["nw"]
end


function read_file_to_vectors(filename::String)
    vectors = []  # Array to hold the vectors
    open(filename, "r") do file
        for line in eachline(file)
            # Split the line by whitespace and convert to numbers
            vector = parse.(Float64, split(line))
            push!(vectors, vector)
        end
    end
    return vectors
end


