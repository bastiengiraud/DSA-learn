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

    ndim += length(bus_gen) # double dimension to account for voltages

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

    optimal_setpoints = []
    original_optimal = []

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
    x_hat = []

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
    #@NLparameter(pm.model, x_hat_p[i = 1:N] == x_hat[i]) # updating the parameter has performance gain vs rebuilding the model
    con_sphere = @constraint(pm.model, sqrt(sum((aux[i]-x_hat[i])^2 for i in 1:N)) <= r) # add constraint on radius

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

            # remove constraint
            delete(pm.model, con_sphere)

            x_hat = sample_ops[j,:] # set x_hat for next iteration
            #JuMP.set_value.(x_hat_p, x_hat)
            # update constraint
            con_sphere = @constraint(pm.model, sqrt(sum((aux[i]-x_hat[i])^2 for i in 1:N)) <= r)

            # compute percentage decrease in volume
            decrease_percentage = (previous_volume - v) / previous_volume

            # Update the previous volume to the new volume
            previous_volume = v

            # Check if the decrease is less than 10%
            if decrease_percentage < stopping_percentage
                counter += 1
                println("Volume decrease of less than $(stopping_percentage*100)%. Current consecutive: ", counter,"/$stopping_iteration", "\n")
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
    in_samples = sample_pol(nb_samples) # draw samples from polytope P
    nb_feasible = 0
    nb_infeasible = 0
    pvpq_feasible = 0
    initial_feasible = 0
    correct_feasible = 0
    tollerance = 1e-4

    pm, N, vars, header = instantiate_system_QCRM(data_tight_tmp, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))
    
    pf_results_feas = []
    load_results_feas = []
    op_info_feas = []

    pf_results_infeas = []
    load_results_infeas = []
    op_info_infeas = []

    # avoid modifying original dataset
    data_opf_verif = deepcopy(data_tight_tmp) 

    sampling_time = @elapsed begin
        for i in 1:nb_samples # loop over number of samples in polytope
            
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

            # dictionary placeholder with OP flags. 1 is feasible, 0 is infeasible
            op_flag = Dict(
                "N0" => 1, # flag for base-case feasibility
                "N0P" => 0.0, # active power violation
                "N0Q" => 0.0, # active power violation
                "N0OV" => 0.0, # amount over voltage violation
                "N0UV" => 0.0, # amount under voltage violation
                "N0L" => 0.0, # amount line flow violation
                "N1" => 1, # flag for N1 feasibility
                "N1OV" => 0.0, # amount over voltage violation
                "N1UV" => 0.0, # amount under voltage violation
                "N1L" => 0.0 # amount line flow violation
            )

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
            pg_vio = 0.0
            qg_vio = 0.0

            for i in 0:(length(multinetwork["nw"])-1)
                # https://github.com/lanl-ansi/PowerModels.jl/blob/30e3d392fd95f1c5a81abd248ac3acc9462211b8/src/prob/pf.jl#L286-L292
                # https://lanl-ansi.github.io/PowerModels.jl/stable/power-flow/#PowerModels.calc_branch_flow_ac -> pf solvers don't produce branch flows
                PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)) # solves non linear pf (unbounded), constraints not enforced
                update_data!(multinetwork["nw"]["$i"], PF_res0["solution"]) # update data with PF results
                flows0 = calc_branch_flow_ac(multinetwork["nw"]["$i"]) # compute branch flows
                update_data!(multinetwork["nw"]["$i"], flows0) # add branch flows
                update_data!(PF_res0["solution"], flows0) # add branch flows to solution
                initial_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res0["solution"], tollerance)

                if i == 0
                    # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
                    if initial_feasibility != true    
                        op_flag["N0"] = 0
                        op_flag["N0P"] += pg_vio
                        op_flag["N0Q"] += qg_vio
                        op_flag["N0OV"] += vm_vio_over
                        op_flag["N0UV"] += vm_vio_under
                        op_flag["N0L"] += sm_vio
                    end
                end

                if i != 0
                    # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
                    if initial_feasibility != true
                        op_flag["N1"] = 0
                        op_flag["N1OV"] += vm_vio_over
                        op_flag["N1UV"] += vm_vio_under
                        op_flag["N1L"] += sm_vio
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
                    "N0P" => 0.0, # active power violation
                    "N0Q" => 0.0, # active power violation
                    "N0OV" => 0.0, # amount over voltage violation
                    "N0UV" => 0.0, # amount under voltage violation
                    "N0L" => 0.0, # amount line flow violation
                    "N1" => 1, # flag for N1 feasibility
                    "N1OV" => 0.0, # amount over voltage violation
                    "N1UV" => 0.0, # amount under voltage violation
                    "N1L" => 0.0 # amount line flow violation
                )

                # solve pvpq
                PF_res1 = nothing
                PVPQ_feasibility = nothing
                vm_vio_over = 0.0
                vm_vio_under = 0.0
                sm_vio = 0.0
                pg_vio = 0.0
                qg_vio = 0.0

                for i in 0:(length(multinetwork["nw"])-1)
                    PF_res1 = adjust_PVPQ(multinetwork["nw"]["$i"], 4)
                    update_data!(multinetwork["nw"]["$i"], PF_res1["solution"]) # update data with PF results
                    flows1 = calc_branch_flow_ac(multinetwork["nw"]["$i"]) # compute branch flows
                    update_data!(multinetwork["nw"]["$i"], flows1) # add branch flows
                    update_data!(PF_res1["solution"], flows1) # add branch flows to solution
                    PVPQ_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res1["solution"], tollerance)

                    if i == 0
                        # if PF_res1["termination_status"] != LOCALLY_SOLVED  || PF_res1["primal_status"] != FEASIBLE_POINT || PF_res1["dual_status"] != FEASIBLE_POINT # PF_res1["termination_status"] == LOCALLY_SOLVED != true && PVPQ_feasibility != true
                        if PVPQ_feasibility != true
                            op_flag["N0"] = 0
                            op_flag["N0P"] += pg_vio
                            op_flag["N0Q"] += qg_vio
                            op_flag["N0OV"] += vm_vio_over
                            op_flag["N0UV"] += vm_vio_under
                            op_flag["N0L"] += sm_vio
                        end
                    end
    
                    if i != 0
                        # if PF_res1["termination_status"] != LOCALLY_SOLVED  || PF_res1["primal_status"] != FEASIBLE_POINT || PF_res1["dual_status"] != FEASIBLE_POINT # PF_res1["termination_status"] == LOCALLY_SOLVED != true && PVPQ_feasibility != true
                        if PVPQ_feasibility != true
                            op_flag["N1"] = 0
                            op_flag["N1OV"] += vm_vio_over
                            op_flag["N1UV"] += vm_vio_under
                            op_flag["N1L"] += sm_vio
                        end
                    end
                end

                if op_flag["N0"] == 1 && op_flag["N1"] == 1 
                    # if pvpq feasible add feasible point
                    nb_feasible += 1
                    pvpq_feasible += 1
                    push!(pf_results_feas, PF_res1["solution"])
                    push!(load_results_feas, data_opf_verif)
                    push!(op_info_feas, op_flag) 
                    println("PV/PQ status:", PF_res1["termination_status"] , "\n")
                    println("PV/PQ feasibility: ", PVPQ_feasibility)
                else
                    # do not add infeasible pvpq
                    # push!(pf_results_infeas, PF_res1["solution"])
                    # push!(load_results_infeas, data_opf_verif)
                    # push!(op_info_infeas, op_flag) 
                    # global nb_infeasible += 1

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
                    @constraint(pm.model, con_sphere, sqrt(sum((aux_new[i]-x_hat[i])^2 for i in 1:N)) <= r) # add constraint on radius

                    @objective(pm.model, Min, r); # the objective is to minimize r

                    # check ac correction feasibility
                    PF_res2 = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
                    acpfcorrect_feasibility = nothing
                    vm_vio_over = 0.0
                    vm_vio_under = 0.0
                    sm_vio = 0.0
                    pg_vio = 0.0
                    qg_vio = 0.0

                    # op placeholder
                    op_flag = Dict(
                        "N0" => 1, # flag for base-case feasibility
                        "N0P" => 0.0, # active power violation
                        "N0Q" => 0.0, # active power violation
                        "N0OV" => 0.0, # amount over voltage violation
                        "N0UV" => 0.0, # amount under voltage violation
                        "N0L" => 0.0, # amount line flow violation
                        "N1" => 1, # flag for N1 feasibility
                        "N1OV" => 0.0, # amount over voltage violation
                        "N1UV" => 0.0, # amount under voltage violation
                        "N1L" => 0.0 # amount line flow violation
                    )


                    for i in 0:(length(multinetwork["nw"])-1)
                        acpfcorrect_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res2["solution"]["nw"]["$i"], tollerance)

                        print(pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio)

                        if i == 0
                            if acpfcorrect_feasibility != true # || PF_res2["termination_status"] == LOCALLY_SOLVED == true  
                                op_flag["N0"] = 0
                                op_flag["N0P"] += pg_vio
                                op_flag["N0Q"] += qg_vio
                                op_flag["N0OV"] += vm_vio_over
                                op_flag["N0UV"] += vm_vio_under
                                op_flag["N0L"] += sm_vio
                            end
                        end
        
                        if i != 0
                            if acpfcorrect_feasibility != true # || PF_res2["termination_status"] == LOCALLY_SOLVED != true  
                                op_flag["N1"] = 0
                                op_flag["N1OV"] += vm_vio_over
                                op_flag["N1UV"] += vm_vio_under
                                op_flag["N1L"] += sm_vio
                            end
                        end
                    end

                    if op_flag["N0"] == 1 && op_flag["N1"] == 1
                        nb_feasible += 1
                        correct_feasible += 1
                        push!(pf_results_feas, PF_res2["solution"]["nw"]["0"])
                        push!(load_results_feas, data_opf_verif)
                        push!(op_info_feas, op_flag)
                        print("acpfcorrect feasibility: ", acpfcorrect_feasibility, "\n")
                        println("ACPF correct status:", PF_res2["termination_status"] , "\n")
                    else
                        # push!(pf_results_infeas, PF_res2["solution"]["nw"]["0"])
                        # push!(load_results_infeas, data_opf_verif)
                        # push!(op_info_infeas, op_flag)
                        # global nb_infeasible += 1
                        print("infeasible acpfcorrect added", "\n")
                        println("PV/PQ status:", PF_res1["termination_status"] )
                        println("ACPF correct status:", PF_res2["termination_status"] )
                    end
                end
            end
    
        end
    end

    print("number of feasible samples: ", nb_feasible, "\n")
    
    # get list of feasible operating points x
    feasible_ops_polytope = []
    for i in 1:length(pf_results_feas[:])
        op = get_Full_OP(data_tight_tmp, pf_results_feas[i], load_results_feas[i])
        push!(feasible_ops_polytope, op)
    end

    # get list of infeasible operating points
    infeasible_ops_polytope = []
    for i in 1:length(pf_results_infeas[:])
        op = get_Full_OP(data_tight_tmp, pf_results_infeas[i], load_results_infeas[i])
        push!(infeasible_ops_polytope, op)
    end

    return feasible_ops_polytope, pf_results_feas, op_info_feas, infeasible_ops_polytope, pf_results_infeas, op_info_infeas, nb_feasible, nb_infeasible, pvpq_feasible, initial_feasible, correct_feasible, sampling_time
end



function sample_mvnd(feasible_ops_polytope, data_tight_tmp, mvnd_samples, contingencies_n1)
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
    
    pf_results_feas = []
    load_results_feas = []
    op_info_feas = []

    pf_results_infeas = []
    load_results_infeas = []
    op_info_infeas = []
    
    nb_feasible_mvnd = 0
    nb_infeasible_mvnd = 0
    initial_feasible_mvnd = 0
    pvpq_feasible_mvnd = 0
    correct_feasible_mvnd = 0
    tollerance = 1e-4

    # avoid modifying original dataset
    data_opf_verif = deepcopy(data_tight_tmp)

    mvnd_sampling_time = @elapsed begin
        for i in eachindex(biased_samples[1,:])
            
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

            print("Current MVND sample number: ", i, "\n")

            # dictionary placeholder with OP flags. 1 is feasible, 0 is infeasible
            op_flag = Dict(
                "N0" => 1, # flag for base-case feasibility
                "N0P" => 0.0, # active power violation
                "N0Q" => 0.0, # active power violation
                "N0OV" => 0.0, # amount over voltage violation
                "N0UV" => 0.0, # amount under voltage violation
                "N0L" => 0.0, # amount line flow violation
                "N1" => 1, # flag for N1 feasibility
                "N1OV" => 0.0, # amount over voltage violation
                "N1UV" => 0.0, # amount under voltage violation
                "N1L" => 0.0 # amount line flow violation
            )

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
            pg_vio = 0.0
            qg_vio = 0.0

            for i in 0:(length(multinetwork["nw"])-1)
                PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
                update_data!(multinetwork["nw"]["$i"], PF_res0["solution"]) # update data with PF results
                flows0 = calc_branch_flow_ac(multinetwork["nw"]["$i"]) # compute branch flows
                update_data!(multinetwork["nw"]["$i"], flows0) # add branch flows
                update_data!(PF_res0["solution"], flows0) # add branch flows to solution
                initial_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res0["solution"], tollerance)

                if i == 0
                    # if PF_res0["termination_status"] != LOCALLY_SOLVED  && PF_res0["primal_status"] != FEASIBLE_POINT && PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
                    if initial_feasibility != true
                        op_flag["N0"] = 0
                        op_flag["N0P"] += pg_vio
                        op_flag["N0Q"] += qg_vio
                        op_flag["N0OV"] += vm_vio_over
                        op_flag["N0UV"] += vm_vio_under
                        op_flag["N0L"] += sm_vio
                    end
                end

                if i != 0
                    # if PF_res0["termination_status"] != LOCALLY_SOLVED  && PF_res0["primal_status"] != FEASIBLE_POINT && PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
                    if initial_feasibility != true
                        op_flag["N1"] = 0
                        op_flag["N1OV"] += vm_vio_over
                        op_flag["N1UV"] += vm_vio_under
                        op_flag["N1L"] += sm_vio
                    end
                end

            end                

            # check feasibility
            if op_flag["N0"] == 1 && op_flag["N1"] == 1
                nb_feasible_mvnd += 1
                initial_feasible_mvnd += 1
                push!(pf_results_feas, PF_res0["solution"])
                push!(load_results_feas, data_opf_verif)
                push!(op_info_feas, op_flag) 
                println("initial status:", PF_res0["termination_status"] , "\n")
            else 
                # add infeasible initial sample
                push!(pf_results_infeas, PF_res0["solution"])
                push!(load_results_infeas, data_opf_verif)
                push!(op_info_infeas, op_flag)
                nb_infeasible_mvnd += 1

                # solve pvpq
                PF_res1 = nothing
                PVPQ_feasibility = nothing
                vm_vio_over = 0.0
                vm_vio_under = 0.0
                sm_vio = 0.0
                pg_vio = 0.0
                qg_vio = 0.0

                # op placeholder
                op_flag = Dict(
                    "N0" => 1, # flag for base-case feasibility
                    "N0P" => 0.0, # active power violation
                    "N0Q" => 0.0, # active power violation
                    "N0OV" => 0.0, # amount over voltage violation
                    "N0UV" => 0.0, # amount under voltage violation
                    "N0L" => 0.0, # amount line flow violation
                    "N1" => 1, # flag for N1 feasibility
                    "N1OV" => 0.0, # amount over voltage violation
                    "N1UV" => 0.0, # amount under voltage violation
                    "N1L" => 0.0 # amount line flow violation
                )

                for i in 0:(length(multinetwork["nw"])-1)
                    PF_res1 = adjust_PVPQ(multinetwork["nw"]["$i"], 4)
                    update_data!(multinetwork["nw"]["$i"], PF_res1["solution"]) # update data with PF results
                    flows1 = calc_branch_flow_ac(multinetwork["nw"]["$i"]) # compute branch flows
                    update_data!(multinetwork["nw"]["$i"], flows1) # add branch flows
                    update_data!(PF_res1["solution"], flows1) # add branch flows to solution
                    PVPQ_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res1["solution"], tollerance)

                    if i == 0
                        # if PF_res1["termination_status"] != LOCALLY_SOLVED  || PF_res1["primal_status"] != FEASIBLE_POINT || PF_res1["dual_status"] != FEASIBLE_POINT # PF_res1["termination_status"] == LOCALLY_SOLVED != true && PVPQ_feasibility != true
                        if PVPQ_feasibility != true
                            op_flag["N0"] = 0
                            op_flag["N0P"] += pg_vio
                            op_flag["N0Q"] += qg_vio
                            op_flag["N0OV"] += vm_vio_over
                            op_flag["N0UV"] += vm_vio_under
                            op_flag["N0L"] += sm_vio
                        end
                    end
    
                    if i != 0
                        #if PF_res1["termination_status"] != LOCALLY_SOLVED  || PF_res1["primal_status"] != FEASIBLE_POINT || PF_res1["dual_status"] != FEASIBLE_POINT # PF_res1["termination_status"] == LOCALLY_SOLVED != true && PVPQ_feasibility != true
                        if PVPQ_feasibility != true
                            op_flag["N1"] = 0
                            op_flag["N1OV"] += vm_vio_over
                            op_flag["N1UV"] += vm_vio_under
                            op_flag["N1L"] += sm_vio
                        end
                    end
                end

                if op_flag["N0"] == 1 && op_flag["N1"] == 1 
                    # if pvpq feasible add feasible point
                    nb_feasible_mvnd += 1
                    pvpq_feasible_mvnd += 1
                    push!(pf_results_feas, PF_res1["solution"])
                    push!(load_results_feas, data_opf_verif)
                    push!(op_info_feas, op_flag) 
                    println("PV/PQ status:", PF_res1["termination_status"] , "\n")
                else
                    # do not add infeasible pvpq
                    # push!(pf_results_infeas, PF_res1["solution"])
                    # push!(load_results_infeas, data_opf_verif)
                    # push!(op_info_infeas, op_flag) 
                    # global nb_infeasible_mvnd += 1

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

                    x_hat = biased_samples[:,i]
                    N = length(vars_new)

                    # additional seperating hyperplanes constraints
                    @variable(pm.model, r); # r is the radius of the hypersphere from the sampled point to the relaxation
                    @variable(pm.model, aux_new[1:N]) # add auxiliary variable (x_opt)
                    @constraint(pm.model, aux_new .== vars_new) # constrain auxiliary variables
                    @constraint(pm.model, con_sphere, sqrt(sum((aux_new[i]-x_hat[i])^2 for i in 1:N)) <= r) # add constraint on radius

                    @objective(pm.model, Min, r); # the objective is to minimize r

                    # check ac correction feasibility
                    PF_res2 = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
                    acpfcorrect_feasibility = nothing
                    vm_vio_over = 0.0
                    vm_vio_under = 0.0
                    sm_vio = 0.0
                    pg_vio = 0.0
                    qg_vio = 0.0

                    # op placeholder
                    op_flag = Dict(
                        "N0" => 1, # flag for base-case feasibility
                        "N0P" => 0.0, # active power violation
                        "N0Q" => 0.0, # active power violation
                        "N0OV" => 0.0, # amount over voltage violation
                        "N0UV" => 0.0, # amount under voltage violation
                        "N0L" => 0.0, # amount line flow violation
                        "N1" => 1, # flag for N1 feasibility
                        "N1OV" => 0.0, # amount over voltage violation
                        "N1UV" => 0.0, # amount under voltage violation
                        "N1L" => 0.0 # amount line flow violation
                    )

                    for i in 0:(length(multinetwork["nw"])-1)
                        acpfcorrect_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res2["solution"]["nw"]["$i"], tollerance)

                        if i == 0
                            if acpfcorrect_feasibility != true #|| PF_res2["termination_status"] == LOCALLY_SOLVED != true  
                                op_flag["N0"] = 0
                                op_flag["N0P"] += pg_vio
                                op_flag["N0Q"] += qg_vio
                                op_flag["N0OV"] += vm_vio_over
                                op_flag["N0UV"] += vm_vio_under
                                op_flag["N0L"] += sm_vio
                            end
                        end
        
                        if i != 0
                            if acpfcorrect_feasibility != true #|| PF_res2["termination_status"] == LOCALLY_SOLVED != true  
                                op_flag["N1"] = 0
                                op_flag["N1OV"] += vm_vio_over
                                op_flag["N1UV"] += vm_vio_under
                                op_flag["N1L"] += sm_vio
                            end
                        end
                    end

                    if op_flag["N0"] == 1 && op_flag["N1"] == 1
                        nb_feasible_mvnd += 1
                        correct_feasible_mvnd += 1
                        push!(pf_results_feas, PF_res2["solution"]["nw"]["0"])
                        push!(load_results_feas, data_opf_verif)
                        push!(op_info_feas, op_flag)
                        print("acpfcorrect feasibility: ", acpfcorrect_feasibility, "\n")
                        println("ACPF correct status:", PF_res2["termination_status"] , "\n")
                    else
                        # push!(pf_results_infeas, PF_res2["solution"]["nw"]["0"])
                        # push!(load_results_infeas, data_opf_verif)
                        # push!(op_info_infeas, op_flag)
                        # global nb_infeasible_mvnd += 1
                        print("infeasible acpfcorrect added", "\n")
                        println("PV/PQ status:", PF_res1["termination_status"] )
                        println("ACPF correct status:", PF_res2["termination_status"] )
                    end
                end
            end
        end


    end

    print("number of feasible mvnd samples: ", nb_feasible_mvnd, "\n")

    feasible_ops_mvnd = []
    for i in 1:length(pf_results_feas[:])
        op = get_Full_OP(data_tight_tmp, pf_results_feas[i], load_results_feas[i])
        push!(feasible_ops_mvnd, op)
    end

    infeasible_ops_mvnd = []
    for i in 1:length(pf_results_infeas[:])
        op = get_Full_OP(data_tight_tmp, pf_results_infeas[i], load_results_infeas[i])
        push!(infeasible_ops_mvnd, op)
    end

    return feasible_ops_mvnd, pf_results_feas, op_info_feas, infeasible_ops_mvnd, pf_results_infeas, op_info_infeas, nb_feasible_mvnd, nb_infeasible_mvnd, pvpq_feasible_mvnd, initial_feasible_mvnd, correct_feasible_mvnd, mvnd_sampling_time

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
                # If contingencies_n1 is empty, inude all branches not in cont_out
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



