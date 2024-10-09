include("acpfcorrect.jl")
include("ssa_module.jl")
include("New_module.jl")
include("obbt_lu.jl")
include("polytope.jl")
include("support.jl")

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



function infeasibility_certif_constraints_check_pgvmpd(network_data , Nb_HP, Nb_insamples, variable_loads)

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
    sample_ops = gen_samples_vectors(Nb_HP, ndim, min_lim, max_lim) # size (Nb_HP, x_hat)

    # create a new optimization model pm with tightened bounds, where you add additional constraints to construct seperating hyperplanes
    pm = instantiate_model(data_tight_tmp, QCRMPowerModel, build_opf_load) # PowerModels.build_opf)

    # initialize input space with the tightened bounds, polytope = P
    create_scaled_pol(ndim, min_lim , max_lim)
    v = comp_vol()

    println("Initial volume : ", v)
    push!(volumes, v)

    # initialize vector of variables
    vars = []
    global x_hat = []

    # get generator variables and add to vars
    slack_gen_idx = get_slack_idx(pm)
    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
        if i != slack_gen_idx && data_tight_tmp["gen"]["$i"]["pmax"] > 0.0
            push!(vars, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
        end
    end

    # uncomment for voltage sampling
    # get voltage magnitude variables and add to vars
    gen_indexes = (unique(map(x -> x["gen_bus"], values(pm.data["gen"]))))
    for g in gen_indexes
        push!(vars, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
    end

    # added
    load_indices = sort(unique(map(x -> x["index"], values(pm.data["load"]))))
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

    HP_time = @elapsed begin
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
                
                println("Volume : ", v)
                push!(volumes, v)
            else
                push!(original_optimal, x_hat)

            end

            global x_hat = sample_ops[j,:] # set x_hat for next iteration
            JuMP.set_value.(x_hat_p, x_hat)

        end
    end

    #--------------------- part 2 ------------------------
    # start sampling from within constructed polytope
    nb_samples = Nb_insamples
    global in_samples = sample_pol(nb_samples) # draw samples from polytope P
    global nb_feasible = 0
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

            PF_res0 = solve_pf(data_opf_verif, ACPPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
            initial_feasibility = check_ac_feasibility(data_tight_tmp, PF_res0, tollerance)
            print("Current sample number: ", i, "\n")
            
            # check feasibility
            if PF_res0["termination_status"] == LOCALLY_SOLVED && initial_feasibility == true 
                global nb_feasible += 1
                global initial_feasible += 1
                push!(pf_results, PF_res0)
                push!(load_results, data_opf_verif)
                println("initial status:", PF_res0["termination_status"] , "\n")
                print("initial feasibility: ", initial_feasibility, "\n")
            else 
                # add infeasible initial sample
                push!(pf_results, PF_res0)
                push!(load_results, data_opf_verif)

                # solve pvpq
                PF_res1 = adjust_PVPQ(data_opf_verif, 8)
                PVPQ_feasibility = check_ac_feasibility(data_tight_tmp, PF_res1, tollerance)

                if PF_res1["termination_status"] == LOCALLY_SOLVED && PVPQ_feasibility == true 
                    # if pvpq feasible add feasible point
                    global nb_feasible += 1
                    global pvpq_feasible += 1
                    push!(pf_results, PF_res1)
                    push!(load_results, data_opf_verif)
                    println("PV/PQ status:", PF_res1["termination_status"] , "\n")
                    print("PVPQ feasibility: ", PVPQ_feasibility, "\n")
                else
                    # add infeasible pvpq
                    push!(pf_results, PF_res1)
                    push!(load_results, data_opf_verif)

                    vars_new = []
                    pm = instantiate_model(data_opf_verif, QCRMPowerModel, build_opf_load)

                    # get generator variables and add to vars
                    slack_gen_idx = get_slack_idx(pm)
                    for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
                        if i != slack_gen_idx && data_tight_tmp["gen"]["$i"]["pmax"] > 0.0
                            push!(vars_new, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
                        end
                    end

                    # uncomment for voltage sampling
                    # get voltage magnitude variables and add to vars
                    gen_indexes = unique(map(x -> x["gen_bus"], values(pm.data["gen"])))
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
                    acpfcorrect_feasibility = check_ac_feasibility(data_tight_tmp, PF_res2, tollerance)

                    if PF_res2["termination_status"] == LOCALLY_SOLVED || acpfcorrect_feasibility == true 
                        global nb_feasible += 1
                        global correct_feasible += 1
                        push!(pf_results, PF_res2)
                        push!(load_results, data_opf_verif)
                        print("acpfcorrect feasibility: ", acpfcorrect_feasibility, "\n")
                        println("ACPF correct status:", PF_res2["termination_status"] , "\n")
                    else
                        push!(pf_results, PF_res2)
                        push!(load_results, data_opf_verif)
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
    for i in 1:length(pf_results) 
        # check for voltage violations 
        vm_vio_over, vm_vio_under = check_vm_violations(data_tight_tmp, pf_results[i]["solution"], tollerance)
        push!(over_array, vm_vio_over) # over voltages
        push!(under_array, vm_vio_under) # under voltages

        # check for generator violations
        pg_vio, qg_vio = check_pg_pq_violations(data_tight_tmp, pf_results[i]["solution"], tollerance)
        push!(pg_array, pg_vio)
        push!(qg_array, qg_vio)

        # check for line violations
        sm_vio = check_flow_violations(data_tight_tmp, pf_results[i]["solution"], tollerance)
        push!(sm_array, sm_vio)

    end

    print("number of feasible samples: ", nb_feasible, "\n")
    
    # count zeros and get index
    _, index_pg = find_zeros_and_indexes(pg_array, tollerance)
    _, index_qg = find_zeros_and_indexes(qg_array, tollerance)
    _, index_sm = find_zeros_and_indexes(sm_array, tollerance)
    _, index_ovi = find_zeros_and_indexes(over_array, tollerance)
    _, index_uvi = find_zeros_and_indexes(under_array, tollerance)

    # if all indices are zero i.e. no violation
    common_elements = [x for x in index_qg if (x in index_sm  && x in index_ovi && x in index_uvi && x in index_pg)]

    # if there is a violation
    rest_of_the_indices = [x for x in eachindex(pf_results) if !(x in common_elements)]

    # get list of feasible operating points x
    OPs = []
    for i in common_elements
        op = get_Full_OP(data_tight_tmp, pf_results[i]["solution"], load_results[i])
        push!(OPs, op)
    end

    # get list of infeasible operating points
    OPs_notFeas = []
    for i in rest_of_the_indices
        op = get_Full_OP(data_tight_tmp, pf_results[i]["solution"], load_results[i])
        push!(OPs_notFeas, op)
    end

    return OPs, pf_results, data_tight_tmp, OPs_notFeas, common_elements, rest_of_the_indices, nb_feasible, pvpq_feasible, initial_feasible, correct_feasible, volumes, HP_time, sampling_time
end



function comment_line_in_file(file_path::String, target_line::String)
    # Read the file
    lines = readlines(file_path)
    
    # Process the lines
    modified_lines = String[]
    for line in lines
        # Check if the line starts with the target string
        if startswith(line, target_line)
            # Comment out the line by adding '%'
            push!(modified_lines, "%" * line)
        else
            # Keep the line as is
            push!(modified_lines, line)
        end
    end
    
    # Write the modified lines back to the file
    open(file_path, "w") do file
        for line in modified_lines
            println(file, line)
        end
    end
end


function get_SetPoint(network_data)
    Pg_SetPoint = Dict{String, Float64}()
    for (i, gen) in network_data["gen"]
        if network_data["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && network_data["gen"][i]["pmax"] > 0.0
            Pg_SetPoint[i] = gen["pg"]
            #push!(Pg_SetPoint, (i, gen["pg"] ))
        end
    end 
    return Pg_SetPoint
end

function update_SetPoint(SetPoint, new_value, nb)
    newSetPoint = deepcopy(SetPoint)
    newSetPoint[nb] = new_value
    return newSetPoint
end


function create_system(network)
    export_matpower("file.m", network)
    file_path = "file.m"
    # The target line to detect and comment out
    target_line = "mpc.multinetwork"
    # Call the function to process the file
    comment_line_in_file(file_path, target_line)
    target_line = "mpc.multiinfrastructure"
    # Call the function to process the file
    comment_line_in_file(file_path, target_line)

    sys = System("file.m") # function from PowerSystems

    return sys
end

function add_dyn_system_basic(sys)
    for g in get_components(Generator, sys)

        #Create the dynamic generator
        case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine_oneDoneQ(),
            shaft = shaft_no_damping(),
            avr = avr_type1(),
            prime_mover = tg_type1(), #tg_none(),
            pss = pss_none(),
            )
        #Attach the dynamic generator to the system
        add_component!(sys, case_gen, g)
    end
    return sys
end


function add_dyn_system_basic_withoutTG(sys)
    for g in get_components(Generator, sys)

        #Create the dynamic generator
        case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine1(),
            shaft = shaft_no_damping(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
            )
        #Attach the dynamic generator to the system
        add_component!(sys, case_gen, g)
    end
    return sys
end



function N_1_step(network_data, OPs_Feas, cont_out)

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

    global flag = false 
    global N_1_feasible = []
    threshold = 10-6

    # get number of variables
    pm, N, vars, header = instantiate_system_QCRM(network_data; variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))

    for i in eachindex(OPs_Feas) # loop over every operating point
        OP = OPs_Feas[i] 
        data_opf_verif = deepcopy(network_data)
        
        for g in eachindex(pg_numbers)
            data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = OP[g] # set generator setpoints to operating point
        end
        for v in eachindex(vm_numbers)
            data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = OP[length(pg_numbers)+v] # set voltages to operating point
        end
        for d in eachindex(pd_numbers)
            data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = OP[length(pg_numbers)+length(vm_numbers)+d] 
            pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
            pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
            sd = pd / pf
            qd = sqrt(sd^2 - pd^2)
            data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
        end

        global network_data = deepcopy(data_opf_verif)

        network_data["area_gens"] = Dict()
        network_data["gen_contingencies"] = []
        network_data["branch_contingencies"] = []
        contingencies = []
        contingencies_gen = []

        for (i, branch) in network_data["branch"]
            if i in cont_out
            continue
            else
                push!(contingencies, (idx=parse(Int,i), label="LINE-$(i)", type="branch"))
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
        
        result = run_c1_scopf_modif(multinetwork, ACPPowerModel, Ipopt.Optimizer)
        for l in 0:(length(multinetwork["nw"])-1)
            update_data!(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
        end

        global result = run_c1_scopf_modif(multinetwork, ACPPowerModel, Ipopt.Optimizer)
        
        

        for l in 0:(length(multinetwork["nw"])-1)
            pg_vio, qg_vio = check_pg_pq_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
            push!(pg_array, pg_vio)
            push!(qg_array, qg_vio)

            sm_vio = check_flow_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
            push!(sm_array, sm_vio)

            vm_over_tmp, vm_under_tmp = check_vm_violations(multinetwork["nw"]["$l"], result["solution"]["nw"]["$l"])
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
            all_below_threshold(pg_array, threshold) &&
            all_below_threshold(qg_array, threshold) &&
            all_below_threshold(sm_array, threshold))
            push!(N_1_feasible, i)
        end
        global over_array = []
        global under_array = []
        global pg_array = []
        global qg_array = []
        global sm_array = []

    end
    return total_sm_array, total_qg_array, total_over_array, total_under_array , total_pg_array, N_1_feasible, result["solution"]["nw"], multinetwork["nw"]
end


function SSS_eval(data_tight_HP, global_OPs, pg_numbers, vm_numbers, pd_numbers)

    global total_damp = []
    global total_dist = []
    global total_eigen = []

    for i in eachindex(global_OPs)   
        data_build = deepcopy(data_tight_HP)
        for g in eachindex(pg_numbers)
            data_build["gen"]["$(pg_numbers[g])"]["pg"] = global_OPs[i][g] 
        end
        for v in eachindex(vm_numbers)
            data_build["bus"]["$(vm_numbers[v])"]["vm"] = global_OPs[i][length(pg_numbers)+v] 
        end
        for d in eachindex(pd_numbers)
            data_build["load"]["$(pd_numbers[d])"]["pd"] = global_OPs[i][length(pg_numbers)+length(vm_numbers)+d] 
            pf = data_build["load"]["$(pd_numbers[d])"]["pf"]
            pd = data_build["load"]["$(pd_numbers[d])"]["pd"]
            sd = pd / pf
            qd = sqrt(sd^2 - pd^2)
            data_build["load"]["$(pd_numbers[d])"]["qd"] = qd
        end

        syst = create_system(data_build)
        add_dyn_system_basic(syst)
        stability = small_signal_module(syst)
        push!(total_damp, stability["damping"])
        push!(total_dist, stability["distance"])
        push!(total_eigen, stability["eigenvalues"])
    
    end
    return total_damp, total_dist, total_eigen
end 


function closest_to_zero_indices(arr, N::Int)
    # Create an array of tuples with absolute value and original index
    abs_with_index = [(abs(val), idx) for (idx, val) in enumerate(arr)]
    
    # Sort the array by absolute value
    sorted_abs_with_index = sort(abs_with_index, by = x -> x[1])
    
    # Get the first N indices from the sorted array
    indices = [sorted_abs_with_index[i][2] for i in 1:N]
    
    return indices
end

function closest_to_zero_indices_no_abs(arr, N::Int)
    # Create an array of tuples with absolute value and original index
    abs_with_index = [(val, idx) for (idx, val) in enumerate(arr)]
    
    # Sort the array by absolute value
    sorted_abs_with_index = sort(abs_with_index, by = x -> x[1])
    
    # Get the first N indices from the sorted array
    indices = [sorted_abs_with_index[i][2] for i in 1:N]
    
    return indices
end




function DW_step(data_tight, index_res, cls_OP, pf_results_prev, distance, alpha)

    file = open("C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/Datasets/DW_file_try.txt", "w")

    global file_nb = 0
    global damping_array = []
    global dist_arr = []
    global plot_damping = []
    global plot_dist = []
    global DW_OPS = []

    gen_keys = collect(keys(data_tight["gen"]))
    gen_slack = get_slack_idx(data_tight)
    filtered_gen_keys = filter(key -> key != gen_slack && data_tight["gen"][key]["pmax"] > 0.0, gen_keys)

    # Get the list of active generators
    active_gens = [data_tight["gen"][key] for key in filtered_gen_keys]

    # Compute the length of active generators
    lgt_active_gen = length(active_gens)

    for i in cls_OP 
        data_build = deepcopy(data_tight)
        update_data!(data_build, pf_results_prev[index_res[i]]["solution"])

        sys_studied = create_system(data_build)
        add_dyn_system_basic(sys_studied)
        stability = small_signal_module(sys_studied)

        println(file, "Initial damping :", stability["damping"])
        println(file, "Initial distance :", stability["distance"])   
        println(file, "Initial eigenvalue :", stability["eigenvalues"])   

        global initial_DR = stability["damping"]
        global damping_array_global = []
        global dist_array_global = []

        
        for DW in 1:15
            # perturb all active generators in positive direction
            for (g, gen) in data_build["gen"]
                if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0 # check if not slack bus and not synchronous condenser
                    network_data_copy = deepcopy(data_build)
                    global updated_value = data_build["gen"]["$g"]["pg"] + 0.2 #* data_build["gen"]["$g"]["pmax"] add perturbation to generator active power

                    if network_data_copy["gen"]["$g"]["pmax"] < updated_value # check if perturbation doesn't violate generator limits
                        network_data_copy["gen"]["$g"]["pg"] = network_data_copy["gen"]["$g"]["pmax"]
                        println(file, "Out of the bounds ")
                    else
                        network_data_copy["gen"]["$g"]["pg"] = updated_value
                    end
                    
                    sys_studied = create_system(network_data_copy)
                    add_dyn_system_basic(sys_studied)
        
                    global file_nb += 1 
                    stability = small_signal_module(sys_studied) # check small signal stability of data with updated generator
                    push!(damping_array, stability["damping"]) # fist half of values in damping_array are all positive perturbations
                    push!(dist_arr, stability["distance"])
        
                end
            end 
        
            # perturb all active generators in negative direction
            for (g, gen) in data_build["gen"]
                if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0
                    network_data_copy = deepcopy(data_build)
                    global updated_value = data_build["gen"]["$g"]["pg"] - 0.2 #* data_build["gen"]["$g"]["pmax"] # perturb generator in other direction

                    if network_data_copy["gen"]["$g"]["pmin"] > updated_value # check if gen limits are not violated
                        network_data_copy["gen"]["$g"]["pg"] = network_data_copy["gen"]["$g"]["pmin"]
                        println(file, "Out of the bounds ")
                    else
                        network_data_copy["gen"]["$g"]["pg"] = updated_value
                    end
                    
                    sys_studied = create_system(network_data_copy)
                    add_dyn_system_basic(sys_studied)
        
                    global file_nb += 1 
                    stability = small_signal_module(sys_studied) # perform small signal stability analysis with negatively perturbed generator
                    push!(damping_array, stability["damping"]) # second half of values in damping_array are all positive perturbations
                    push!(dist_arr, stability["distance"])
        
                end
            end 

            index_gen_disturb = argmin(abs.(damping_array)) # return index of smallest value i.e. smallest damping coefficient
            sorted_index = sortperm(abs.(damping_array)) # sort the indices
            min_damp = minimum(damping_array) # return the smallest damping coefficient
            
            # define step size of directed walks based on distance
            if abs(min_damp) > distance[1]
                step_size = alpha[1]
            elseif distance[2] < abs(min_damp)  < distance[1]
                step_size = alpha[2]
            elseif distance[3] < abs(min_damp) < distance[2]
                step_size = alpha[3] #0.15
            elseif abs(min_damp) < distance[3]
                step_size = alpha[4] #0.15    
            end

            # determine the direction of the directed walk -> due to positive or negative generator perturbation
            println(file, "Directed walk number: $DW __________________")
            sign_DW = 1
            if index_gen_disturb > (lgt_active_gen) # if the min damping is in the second half of the array, the min damping is due to a negative perturbation
                sign_DW = -1
                index_gen_disturb -= lgt_active_gen # subtract number of generators to obtain generator index
            end    
            
            # check of the current generator setpoint violates constraints, if so, take the next damping value
            if data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] <= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmin"] || data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] >= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmax"]
                index_gen_disturb = sorted_index[2]
                if index_gen_disturb > (lgt_active_gen)
                    sign_DW = -1
                    index_gen_disturb -= lgt_active_gen
                end    

                # do same check for the next value
                if data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] <= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmin"] || data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pg"] >= data_build["gen"][filtered_gen_keys[index_gen_disturb]]["pmax"]
                    index_gen_disturb = sorted_index[3]
                    if index_gen_disturb > (lgt_active_gen)
                        sign_DW = -1
                        index_gen_disturb -= lgt_active_gen
                    end    
                    
                end

            end

            # start directed walks
            sP = get_SetPoint(data_build) # get the generator setpoints
            println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
            index_gen_disturb = filtered_gen_keys[index_gen_disturb]
            newOP_PG = data_build["gen"][index_gen_disturb]["pg"] + sign_DW * 0.2 #* data_build["gen"][index_gen_disturb]["pmax"] # update the generator setpoints
            new_SP = update_SetPoint(sP, newOP_PG, index_gen_disturb) # update the new setpoints
            println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
            
            grad = (collect(values(new_SP)) - collect(values(sP)))/norm(collect(values(new_SP)) - collect(values(sP)))
            println(file, "grad: ", join(grad, ", "))

            delta_gen =  step_size*data_build["gen"][index_gen_disturb]["pg"].*(grad) #.*collect(values(new_SP))) # step generator
            println(file, "delta generator: ", join(delta_gen, ", "))
            # println(file, data_build["gen"] )

            update_gen_network(data_build, delta_gen) # update the generator setpoint i.e. Pgn+1 = Pgn - delta_gen
            println(file, "new generator setpoint: ", join(get_SetPoint(data_build), ", "))

            # check small signal stability of perturbered generator and obtain stability indices
            sys_studied = create_system(data_build)
            add_dyn_system_basic(sys_studied)
            stability = small_signal_module(sys_studied)

            println(file, "damping of new setpoint: ", stability["damping"])
            println(file, "distance of new setpoint: ", stability["distance"])
            println(file, "eigenvalues of new setpoint: ", stability["eigenvalues"])

            push!(plot_damping, stability["damping"])
            push!(plot_dist, stability["distance"])
            push!(damping_array_global , damping_array)

            global damping_array = []
            push!(dist_array_global , dist_arr)

            global dist_arr = []
            OP_tmp = get_Full_OP(data_build, data_build, data_build)
            push!(DW_OPS, OP_tmp)
            push!(DW_OPS, (stability["damping"], stability["distance"]))


            if initial_DR < stability["damping"] # stop directed walks if 
                break
            end
        end

        println(file, dist_array_global)
        println(file, damping_array_global)
        dist_array_global = []
        damping_array_global = []
        
    end
    close(file)

    return DW_OPS 

end


function sample_MVND(OP_HP, network_basic, data_tight, Nb_ops)
    mean_vec = mean(OP_HP, dims=1)

    cov_matrix = cov(OP_HP)
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
    biased_samples = rand(mvn, Nb_ops)
    pm, N, vars, header = instantiate_system_QCRM(data_tight, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))
    pf_results = []
    load_results = []
    
    global nb_feasible_mvnd = 0
    tollerance = 1e-9

    for i in eachindex(biased_samples[1,:])
        data_opf_verif = deepcopy(data_tight)
        
        for g in eachindex(pg_numbers)
            data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = biased_samples[g,i] #optimal_setpoints[i][g] #in_samples[g,i]
        end
        for v in eachindex(vm_numbers)
            data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = biased_samples[length(pg_numbers)+v,i] #optimal_setpoints[i][length(pg_numbers)+v] #in_samples[length(pg_numbers)+v,i]
        end
        for d in eachindex(pd_numbers)
            data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = biased_samples[length(pg_numbers)+length(vm_numbers)+d,i] 
            pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
            pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
            sd = pd / pf
            qd = sqrt(sd^2 - pd^2)
            data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
        end

        # solve pf, and if reactive power violation, change to PV bus and solve pf again
        PF_res1 = adjust_PVPQ(data_opf_verif, 8)
        PVPQ_feasibility = check_ac_feasibility(data_tight, PF_res1, tollerance)

        # check feasibility
        if PF_res1["termination_status"] == LOCALLY_SOLVED && PVPQ_feasibility == true 
            global nb_feasible_mvnd += 1
            push!(pf_results, PF_res1)
            push!(load_results, data_opf_verif)
            println("PV/PQ status:", PF_res1["termination_status"] , "\n")
            print("PVPQ feasibility: ", PVPQ_feasibility, "\n")
        else
            # additional acpf correction! 
            # push pg_cmd to data
            push_generator_setpoints(data_opf_verif)

            # solve correction
            solver = optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
            m = Model()
            set_silent(m)
            PF_res2 = solve_model(data_opf_verif, ACPPowerModel, solver, build_acpf_correction ,jump_model=m) 
            acpfcorrect_feasibility = check_ac_feasibility(data_tight, PF_res2, tollerance)

            if PF_res2["termination_status"] == LOCALLY_SOLVED || acpfcorrect_feasibility == true 
                global nb_feasible_mvnd += 1
                push!(pf_results, PF_res2)
                push!(load_results, data_opf_verif)
                print("acpfcorrect feasibility: ", acpfcorrect_feasibility, "\n")
                println("ACPF correct status:", PF_res2["termination_status"] , "\n")
            else
                push!(pf_results, PF_res2)
                push!(load_results, data_opf_verif)
                print("infeasible acpfcorrect added", "\n")
                println("PV/PQ status:", PF_res1["termination_status"] )
                println("ACPF correct status:", PF_res2["termination_status"] )
            end

        end


    end

    global over_array = []
    global under_array = []
    global pg_array = []
    global qg_array = []
    global sm_array = []
    global res_arr = []

    # check for feasibility
    for i in 1:length(pf_results)
        vm_vio_over, vm_vio_under = check_vm_violations(network_basic, pf_results[i]["solution"], tollerance)
        push!(over_array, vm_vio_over)
        push!(under_array, vm_vio_under)

        pg_vio, qg_vio = check_pg_pq_violations(network_basic, pf_results[i]["solution"], tollerance)
        push!(pg_array, pg_vio)
        push!(qg_array, qg_vio)

        sm_vio = check_flow_violations(network_basic, pf_results[i]["solution"], tollerance)
        push!(sm_array, sm_vio)

    end

    nb_pg, index_pg = find_zeros_and_indexes(pg_array, tollerance)
    nb_qg, index_qg = find_zeros_and_indexes(qg_array, tollerance)
    nb_sm, index_sm = find_zeros_and_indexes(sm_array, tollerance)
    nb_ovi, index_ovi = find_zeros_and_indexes(over_array, tollerance)
    nb_uvi, index_uvi = find_zeros_and_indexes(under_array, tollerance)

    print("number of feasible mvnd samples: ", nb_feasible_mvnd, "\n")

    common_elements = [x for x in index_qg if (x in index_sm  && x in index_ovi && x in index_uvi && x in index_pg)]
    rest_of_the_indices = [x for x in eachindex(pf_results) if !(x in common_elements)]

    OPs = []
    for i in common_elements
        op = get_Full_OP(data_tight, pf_results[i]["solution"], load_results[i])
        push!(OPs, op)
    end

    OPs_notFeas = []
    for i in rest_of_the_indices
        op = get_Full_OP(data_tight, pf_results[i]["solution"], load_results[i])
        push!(OPs_notFeas, op)
    end

    return OPs, OPs_notFeas, nb_feasible_mvnd

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


