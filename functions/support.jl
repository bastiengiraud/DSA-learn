


function header_full_op(data_model, header, pg_numbers, vm_numbers)
    op_header = deepcopy(header)
    gen_vars = length(pg_numbers)
    volt_vars = length(vm_numbers)
    op_header = op_header[1][1:(gen_vars+volt_vars)]
    loads = length(data_model["load"])

    for i in 1:loads
        push!(op_header, "PD$i")
    end

    return [op_header]
end

function get_slack_idx(power_model::AbstractPowerModel)
    # Slack type == 3
    bus_idx =  [k for (k,v) in power_model.data["bus"] if v["bus_type"] == 3]
    bus_idx = parse(Int64,bus_idx[1])
    gen_idx = [k for (k,v) in power_model.data["gen"] if v["gen_bus"] == bus_idx]
    return parse(Int64, gen_idx[1])
end

function get_slack_idx_mn(power_model::AbstractPowerModel)
    # Slack type == 3
    bus_idx =  [k for (k,v) in power_model.data["nw"]["0"]["bus"] if v["bus_type"] == 3]
    bus_idx = parse(Int64,bus_idx[1])
    gen_idx = [k for (k,v) in power_model.data["nw"]["0"]["gen"] if v["gen_bus"] == bus_idx]
    return parse(Int64, gen_idx[1])
end

function get_slack_idx(network)
    # Slack type == 3
    bus_idx =  [k for (k,v) in network["bus"] if v["bus_type"] == 3]
    bus_idx = parse(Int64,bus_idx[1])
    gen_idx = [k for (k,v) in network["gen"] if v["gen_bus"] == bus_idx]
    return  gen_idx[1]
end


function get_header_names(variables::AbstractArray)
    header = [] # used as header in x_opt results csv
    push!(header, JuMP.name.(variables))
    header = map((x) -> uppercase(replace(x, r"0_|\[|\]" => "")),header[1]) # prettify
    return [header]::AbstractArray
end

function gen_samples(n_samples,n_dimensions,level_min, level_max)
    scaling_list = [(level_min, level_max) for _ in 1:n_dimensions]
    plan, _ = LHCoptim(n_samples,n_dimensions,1000)
    scaled_plan = scaleLHC(plan,scaling_list)
    return scaled_plan
end

function gen_samples_vectors(n_samples, n_dimensions, level_min, level_max)
    scaling_list = [(level_min[i], level_max[i]) for i in eachindex(level_min)]
    plan, _ = LHCoptim(n_samples, n_dimensions, 1000)
    scaled_plan = scaleLHC(plan, scaling_list)
    return scaled_plan
end

function update_all_demand(network_data, sample_array)
    for i in 1:length(sample_array)
        network_data["load"]["$i"]["pd"] *= sample_array[i]
        network_data["load"]["$i"]["qd"] *= sample_array[i]
    end
    return network_data
end

function scale_all_demand(network_data, scaling)
    nb_load = length(network_data["load"])
    for i in 1:nb_load
        network_data["load"]["$i"]["pd"] *= scaling
        network_data["load"]["$i"]["qd"] *= scaling
    end
    return network_data
end

function update_all_demand_reactive(network_data, sample_array)
    for i in 1:length(sample_array)
        network_data["load"]["$i"]["qd"] *= sample_array[i]
    end
    return network_data
end


function update_all_limits_reactive(network_data, sample_array)
    for i in 1:length(sample_array)
        network_data["gen"]["$i"]["qmax"] *= sample_array[i]
        network_data["gen"]["$i"]["qmin"] *= sample_array[i]
    end
    return network_data
end

function OPF_feasible_samples(sample_array, network_data)
    global feasible = 0
    global infeasible = 0 
    global non_feasible_index = []
    global feasible_index = []
    for i in 1:length(sample_array[:,1])
        network_data_tmp = deepcopy(network_data)
        demand_profile = sample_array[i,:]
        update_all_demand(network_data_tmp, demand_profile)
        results_QCOPF = solve_opf(network_data_tmp,QCLSPowerModel,optimizer_with_attributes(solver, "print_level" => 0))
        # Check feasibility
        feasible_iteration = 0
        if results_QCOPF["termination_status"] == LOCALLY_SOLVED
            feasible_iteration = 1 
            push!(feasible_index,i)
        else
            println(results_QCOPF["termination_status"])
            push!(non_feasible_index,i)
        end

        global feasible += feasible_iteration
        global infeasible += (1 - feasible_iteration)
    end
    return feasible, infeasible, non_feasible_index, feasible_index
end

function OPF_feasible_samples_reactive(sample_array, network_data)
    global feasible = 0
    global infeasible = 0 
    global non_feasible_index = []
    global feasible_index = []
    for i in 1:length(sample_array[:,1])
        network_data_tmp = deepcopy(network_data)
        demand_profile = sample_array[i,:]
        update_all_demand_reactive(network_data_tmp, demand_profile)
        results_QCOPF = solve_opf(network_data_tmp,QCLSPowerModel,optimizer_with_attributes(solver, "print_level" => 0))
        # Check feasibility
        feasible_iteration = 0
        if results_QCOPF["termination_status"] == LOCALLY_SOLVED
            feasible_iteration = 1 
            push!(feasible_index,i)
        else
            println(results_QCOPF["termination_status"])
            push!(non_feasible_index,i)
        end

        global feasible += feasible_iteration
        global infeasible += (1 - feasible_iteration)
    end
    return feasible, infeasible, non_feasible_index, feasible_index
end


function extract_number_and_type(variable_names::Vector{String})
    pg_numbers = Int[]
    vm_numbers = Int[]
    pd_numbers = Int[]
    for variable_name in variable_names
        match_result = match(r"(PG|VM|PD)(\d+)", variable_name)
        if match_result !== nothing
            variable_type = match_result.captures[1]
            variable_number = parse(Int, match_result.captures[2])
            if variable_type == "PG"
                push!(pg_numbers, variable_number)
            elseif variable_type == "VM"
                push!(vm_numbers, variable_number)
            elseif variable_type == "PD"
                push!(pd_numbers, variable_number)
            end
        end
    end
    return pg_numbers, vm_numbers, pd_numbers
end

function extract_number_and_type(variable_name::String)
    match_result = match(r"(PG|VM|PD)(\d+)", variable_name)
    if match_result !== nothing
        variable_type = match_result.captures[1]
        variable_number = parse(Int, match_result.captures[2])
        return variable_type, variable_number
    else
        return nothing, nothing
    end
end


function randomize_values(n::Int, values::Vector{Float64})
    if isempty(values)
        throw(ArgumentError("Input vector cannot be empty"))
    end
    
    min_val = minimum(values)
    max_val = maximum(values)
    
    random_values = [rand(length(values)) .* (max_val - min_val) .+ min_val for _ in 1:n]
    
    return random_values
end


function update_SetPoint(SetPoint, new_value, nb)
    newSetPoint = deepcopy(SetPoint)
    newSetPoint[nb] = new_value
    return newSetPoint
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



function update_gen_network(network, new_gen)
    ind = 1
    for (i, gen) in network["gen"]
        if network["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && network["gen"][i]["pmax"] > 0.0
            gen["pg"] += new_gen[ind]
            if gen["pg"] < gen["pmin"]
                gen["pg"] = gen["pmin"]
            elseif gen["pg"] > gen["pmax"]
                gen["pg"] = gen["pmax"]
            end

            ind += 1
        end
    end 
end

function update_bus_network(network, active_gens_buses ,new_bus)
    ind = 1
    for i in active_gens_buses
        network["bus"]["$i"]["vm"] += new_bus[ind]
        if network["bus"]["$i"]["vm"] < network["bus"]["$i"]["vmin"]
            network["bus"]["$i"]["vm"] = network["bus"]["$i"]["vmin"]
        elseif network["bus"]["$i"]["vm"] > network["bus"]["$i"]["vmax"]
            network["bus"]["$i"]["vm"] = network["bus"]["$i"]["vmax"]
        end

        ind += 1
    
    end 
end

function get_SetPoint(network_data,active_gens_buses)
    Pg_SetPoint = Dict{String, Float64}()
    for (i, gen) in network_data["gen"]
        if network_data["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && network_data["gen"][i]["pmax"] > 0.0
            Pg_SetPoint[i] = gen["pg"]
            #push!(Pg_SetPoint, (i, gen["pg"] ))
        end
    end 
    return Pg_SetPoint
end

function get_SetPoint_V(network_data, active_gens_buses)
    Pg_SetPoint = Dict{Int, Float64}()
    for (bus_key, bus_data) in network_data["bus"]
        bus_idx = parse(Int, bus_key)
        if bus_idx in active_gens_buses
            Pg_SetPoint[bus_idx] = bus_data["vm"]
        end
    end
    return Pg_SetPoint
end



function get_OP_from_System(network)
    slack_gen_idx = get_slack_idx(network)
    OP_x = []
    for (g, gen) in network["gen"]
        if g != slack_gen_idx && gen["pmax"] > 0.0
            push!(OP_x, gen["pg"] )
        end
    end
    gen_indexes = unique(map(x -> x["gen_bus"], values(network["gen"])))
    for b in gen_indexes
        push!(OP_x, network["bus"]["$b"]["vm"])
    end
    return OP_x
end

function all_below_threshold(arr, threshold)
    return all(x -> x < threshold, arr)
end

function get_Full_OP(network_data_basic, network_data, load_results)
    OP = []
    slack_gen_idx = get_slack_idx(network_data_basic)

    for i=1:length(network_data["gen"])
        if i != parse(Int,slack_gen_idx) && network_data_basic["gen"]["$i"]["pmax"] > 0.0
            push!(OP, network_data["gen"]["$i"]["pg"])
        elseif i == slack_gen_idx
            continue
        end
    end

    gen_indexes = unique(x["gen_bus"] for x in values(network_data_basic["gen"]))
    for g in gen_indexes
        push!(OP, network_data["bus"]["$g"]["vm"])
    end

    load_indexes = unique(x["load_bus"] for x in values(load_results["load"]))
    for d in 1:length(load_indexes)
        push!(OP, load_results["load"]["$d"]["pd"])
    end

    return OP
end 


function find_zeros_and_indexes(vector, tollerance)
    # Identify the indexes where the value is within the range [-threshold, threshold]
    zero_indexes = findall(x -> abs(x) <= tollerance, vector)
    # zero_indexes = findall(x -> x == 0, vector)
    
    # Count the number of zeros (or near zeros)
    zero_count = length(zero_indexes)
    
    return zero_count, zero_indexes
end


function check_vm_violations(network, solution, tollerance; vm_digits = 6)
    vm_vio_over = 0.0
    vm_vio_under = 0.0
    nb_vm = length(network["bus"])

        for (i,bus) in network["bus"]
            if bus["bus_type"] != 4
                bus_sol = solution["bus"][i]

                # helps to account for minor errors in equality constraints
                sol_val = round(bus_sol["vm"], digits=vm_digits)

                if sol_val < bus["vmin"] - tollerance
                    vm_vio_under += bus["vmin"] - sol_val

                end
                if sol_val > bus["vmax"] + tollerance
                    vm_vio_over += sol_val - bus["vmax"]

                end

            end
        end
        vm_vio_over = vm_vio_over#/nb_vm
        vm_vio_under = vm_vio_under#/nb_vm
        return vm_vio_over, vm_vio_under 
end


function check_vm_violations_mn(network, solution, nw, tollerance; vm_digits = 6)
    vm_vio_over = 0.0
    vm_vio_under = 0.0
    nb_vm = length(network["bus"])

        for (i,bus) in network["bus"]
            if bus["bus_type"] != 4
                bus_sol = solution["nw"]["$nw"]["bus"][i]

                # helps to account for minor errors in equality constraints
                sol_val = round(bus_sol["vm"], digits=vm_digits)

                if sol_val < bus["vmin"] - tollerance
                    vm_vio_under += bus["vmin"] - sol_val

                end
                if sol_val > bus["vmax"] + tollerance
                    vm_vio_over += sol_val - bus["vmax"]

                end

            end
        end
        vm_vio_over = vm_vio_over#/nb_vm
        vm_vio_under = vm_vio_under#/nb_vm
        return vm_vio_over, vm_vio_under 
end

function check_vm_violations_bis(network, solution; vm_digits = 6)
    vm_vio_over = 0.0
    vm_vio_under = 0.0
    nb_vm = length(network["bus"])

    for (i,bus) in network["bus"]
        if bus["bus_type"] != 4
            bus_sol = solution["bus"][i]

            # helps to account for minor errors in equality constraints
            sol_val = round(bus_sol["vm"], digits=vm_digits)

            if sol_val < bus["vmin"]
                vm_vio_under += bus["vmin"] - sol_val
            end
            if sol_val > bus["vmax"]
                vm_vio_over += sol_val - bus["vmax"]
            end

        end
    end
    vm_vio_over = vm_vio_over
    vm_vio_under = vm_vio_under
    return vm_vio_over, vm_vio_under 
end


function check_pg_pq_violations(network, solution, tollerance)
    pg_vio = 0.0
    qg_vio = 0.0
    nb_gen = length(network["gen"])
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0 && network["bus"]["$(gen["gen_bus"])"]["bus_type"] != 4
            gen_sol = solution["gen"][i]

            if gen_sol["pg"] < gen["pmin"] - tollerance
                pg_vio += abs(gen["pmin"] - gen_sol["pg"])
                # print("violation at gen: ", i, " Pg is: ", gen_sol["pg"], "\n")
                # if gen["pmax"] > 0
                #     pg_vio += ((gen["pmin"] - gen_sol["pg"])/gen["pmax"])*100
                # else
                #     pg_vio += (gen["pmin"] - gen_sol["pg"])
                # end
            end
            if gen_sol["pg"] > gen["pmax"] + tollerance
                pg_vio += abs(gen_sol["pg"] - gen["pmax"])
                # print("violation at gen: ", i, " Pg is: ", gen_sol["pg"], "\n")
                # if gen["pmax"] > 0
                #     pg_vio += ((gen_sol["pg"] - gen["pmax"])/gen["pmax"])*100
                # else
                #     pg_vio += (gen["pmin"] - gen_sol["pg"])
                # end
            end

            if gen_sol["qg"] < gen["qmin"] - tollerance
                qg_vio += (abs(gen["qmin"] - gen_sol["qg"]))
                #qg_vio += (abs(gen["qmin"] - gen_sol["qg"])/gen["qmax"])*100
                # print("violation at gen: ", i, " Qg is: ", gen_sol["qg"], "\n")
            end
            if gen_sol["qg"] > gen["qmax"] + tollerance
                qg_vio += (abs(gen_sol["qg"] - gen["qmax"]))
                # qg_vio += ((gen_sol["qg"] - gen["qmax"])/gen["qmax"])*100
                # print("violation at gen: ", i, " Qg is: ", gen_sol["qg"], "\n")
            end
        end
    end
    pg_vio = pg_vio#/nb_gen
    qg_vio = qg_vio#/nb_gen
    return pg_vio, qg_vio

end


function check_pg_pq_violations_mn(network, solution, nw, tollerance)
    pg_vio = 0.0
    qg_vio = 0.0
    nb_gen = length(network["gen"])
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0 && network["bus"]["$(gen["gen_bus"])"]["bus_type"] != 4
            gen_sol = solution["nw"]["$nw"]["gen"][i]

            if gen_sol["pg"] < gen["pmin"] - tollerance
                if gen["pmax"] > 0
                    pg_vio += ((gen["pmin"] - gen_sol["pg"])/gen["pmax"])*100
                else
                    pg_vio += (gen["pmin"] - gen_sol["pg"])
                end
            end
            if gen_sol["pg"] > gen["pmax"] + tollerance
                if gen["pmax"] > 0
                    pg_vio += ((gen_sol["pg"] - gen["pmax"])/gen["pmax"])*100
                else
                    pg_vio += (gen["pmin"] - gen_sol["pg"])
                end
            end

            if gen_sol["qg"] < gen["qmin"] - tollerance
                qg_vio += (abs(gen["qmin"] - gen_sol["qg"])/gen["qmax"])*100
            end
            if gen_sol["qg"] > gen["qmax"] + tollerance
                qg_vio += ((gen_sol["qg"] - gen["qmax"])/gen["qmax"])*100
            end
        end
    end
    pg_vio = pg_vio#/nb_gen
    qg_vio = qg_vio#/nb_gen
    return pg_vio, qg_vio

end

function check_flow_violations(network, solution, tollerance)
    nb_branch = length(network["branch"])
    sm_vio = 0.0 #NaN
    if haskey(solution, "branch")
        sm_vio = 0.0
        for (i,branch) in network["branch"]
            if branch["br_status"] != 0

                if !haskey(solution["branch"], i)
                    continue  # if there is no key (contingency), skip key
                end
                
                branch_sol = solution["branch"][i]

                s_fr = abs(branch_sol["pf"])
                s_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                rating = branch["rate_a"]

                if s_fr > rating + tollerance
                    sm_vio += abs(s_fr - rating) #((s_fr - rating)/rating)*100
                end
                if s_to > rating + tollerance
                    sm_vio += abs(s_to - rating) #((s_to - rating)/rating)*100
                end
            end
        end

    else
        throw(ErrorException("There are no branch flows in the solution. Can't check AC feasiblity."))

    end

    sm_vio = sm_vio#/nb_branch

    return sm_vio
end

function check_flow_violations_mn(network, solution, nw, tollerance)
    nb_branch = length(network["branch"])
    sm_vio = 0.0 #NaN
    if haskey(solution, "branch")
        sm_vio = 0.0
        for (i,branch) in network["branch"]
            if branch["br_status"] != 0
                branch_sol = solution["nw"]["$nw"]["branch"][i]

                s_fr = abs(branch_sol["pf"])
                s_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                rating = branch["rate_a"]

                if s_fr > rating + tollerance
                    sm_vio += abs(s_fr - rating) #((s_fr - rating)/rating)*100
                end
                if s_to > rating + tollerance
                    sm_vio += abs(s_to - rating) #((s_to - rating)/rating)*100
                end
            end
        end
    end

    sm_vio = sm_vio#/nb_branch

    return sm_vio
end


function check_flow_violations_bis(network, solution)

    if haskey(network, "branch")
        sm_vio = []
        s_fr_arr = []
        s_to_arr = []
        for (i,branch) in network["branch"]
            if isnan(branch["qf"]) &&  isnan(branch["qf"]) && isnan(branch["qt"]) && isnan(branch["qt"])
                continue 
            else
                s_fr = abs(branch["pf"])
                s_to = abs(branch["pt"])
                
                if !isnan(branch["qf"]) && !isnan(branch["qt"])
                    s_fr = sqrt(branch["pf"]^2 + branch["qf"]^2)
                    s_to = sqrt(branch["pt"]^2 + branch["qt"]^2)
                end
                push!(s_fr_arr, s_fr)
                push!(s_to_arr, s_to)

                rating = network["branch"]["$i"]["rate_c"]

                if s_fr > rating
                    push!(sm_vio,i)
                    push!(sm_vio, ((s_fr - rating)/rating)*100)
                end
                if s_to > rating
                    push!(sm_vio,i)
                    push!(sm_vio, ((s_to - rating)/rating)*100)
                end

            end
        end
    end
    sm_vio = sm_vio/nb_branch
    return sm_vio, s_fr_arr, s_to_arr
end


function check_ac_feasibility(data_tight_tmp, result, tollerance)
    # check for voltage violations
    vm_vio_over, vm_vio_under = check_vm_violations(data_tight_tmp, result, tollerance)
    
    # check for generator violations
    pg_vio, qg_vio = check_pg_pq_violations(data_tight_tmp, result, tollerance)
    
    # check for line violations
    sm_vio = check_flow_violations(data_tight_tmp, result, tollerance)

    if (vm_vio_over + vm_vio_under + pg_vio + qg_vio + sm_vio) == 0
        return true, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio
    else 
        return false, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio
    end

end

function check_ac_feasibility_mn(data_tight_tmp, result, tollerance)
    nb_contingencies = length(result["solution"]["nw"])
    feasible = 0

    for nw in 0:(nb_contingencies-1)
        # check for voltage violations
        vm_vio_over, vm_vio_under = check_vm_violations_mn(data_tight_tmp, result["solution"], nw, tollerance)
        
        # check for generator violations
        pg_vio, qg_vio = check_pg_pq_violations_mn(data_tight_tmp, result["solution"], nw, tollerance)
        
        # check for line violations
        sm_vio = check_flow_violations_mn(data_tight_tmp, result["solution"], nw, tollerance)

        if (vm_vio_over + vm_vio_under + pg_vio + qg_vio + sm_vio) == 0
            feasible += 0
        else 
            feasible += 1
        end
    
    end

    if feasible == 0
        return true
    else
        return false
    end

end


function where_ac_violation(data_tight_tmp, result, tollerance)
    # check for voltage violations
    vm_vio_over, vm_vio_under = check_vm_violations(data_tight_tmp, result["solution"], tollerance)
    
    # check for generator violations
    pg_vio, qg_vio = check_pg_pq_violations(data_tight_tmp, result["solution"], tollerance)
    
    # check for line violations
    sm_vio = check_flow_violations(data_tight_tmp, result["solution"], tollerance)

    return vm_vio_over, vm_vio_under, pg_vio, qg_vio, sm_vio

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



function update_active_reactive_power_data!(network::Dict{String,<:Any}, data::Dict{String,<:Any}; branch_flow=false)
    for (i,bus) in data["bus"]
        nw_bus = network["bus"][i]
        nw_bus["va"] = bus["va"]
    end

    for (i,gen) in data["gen"]
        nw_gen = network["gen"][i]
        nw_gen["pg"] = gen["pg"]
        nw_gen["qg"] = gen["qg"]
    end

    if branch_flow
        for (i,branch) in data["branch"]
            nw_branch = network["branch"][i]
            nw_branch["pf"] = branch["pf"]
            nw_branch["pt"] = branch["pt"]
        end
    end
end

#using Filesystem
function clear_temp_folder(temp_folder::String)
    for entry in readdir(temp_folder)
        file_path = joinpath(temp_folder, entry)
        try
            isdir(file_path) ? rm(file_path; recursive=true) : rm(file_path)
        catch e
            @warn "Failed to remove $file_path: $e"
        end
    end
end


function check_initialization()
    # call error if DWs is true but sss is not true
    if Initialize.sss_analysis != true && Initialize.directed_walks == true
        throw(ErrorException("You can't perform directed walks without small-signal stability analysis."))
    end
end


function check_initialization_lhc()
    # call error if DWs is true or mvnd is true
    if Initialize.mvnd_sampling == true || Initialize.directed_walks == true
        throw(ErrorException("The LHC baseline doesn't do importance sampling or directed walks."))
    end
end


function check_initialization_imp()
    # call error if DWs is true or mvnd is false
    if Initialize.directed_walks == true
        throw(ErrorException("The Importance sampling baseline doesn't do directed walks."))
    elseif Initialize.mvnd_sampling != true 
        throw(ErrorException("The Importance sampling baseline needs mvnd sampling."))
    end
end


function clear_custom_temp_folder()
    temp_folder = ENV["JULIA_TEMP"]
    
    # Check if the folder exists
    if isdir(temp_folder)
        # Remove all files in the folder
        for file in readdir(temp_folder)
            rm(joinpath(temp_folder, file), recursive=true)
        end
        println("Custom temp folder cleared.")
    else
        println("Custom temp folder does not exist.")
    end
end

function clean_temp_files()
    temp_directory = tempdir()  # Get the system's temporary directory
    for file in readdir(temp_directory)  # List all files in the temp directory
        if startswith(file, "jl_")  # Check if the file starts with 'jl_'
            try
                rm(joinpath(temp_directory, file))  # Remove the file
            catch e
                @warn "Could not remove file $file: $e"
            end
        end
    end
end


function clean_full_temp_folder()
    temp_dir = tempdir()  # Get the system's temporary directory
    println("Cleaning temp folder: $temp_dir")
    
    # List all files and directories in the temp folder
    for entry in readdir(temp_dir, join=true)
        try
            if isdir(entry)
                rm(entry, recursive=true)  # Remove directories and their contents
            else
                rm(entry)  # Remove files
            end
        catch e
            @warn "Failed to remove $entry" exception=(e, catch_backtrace())
        end
    end
    println("Temp folder cleaned successfully!")
end
