
function constraint_gen_setpoint_reactive_max(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)
    constraint_gen_setpoint_reactive(pm, nw, gen["index"], gen["qmax"]) # qg[i] = qmax
end

function constraint_gen_setpoint_reactive_min(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    gen = ref(pm, nw, :gen, i)
    constraint_gen_setpoint_reactive(pm, nw, gen["index"], gen["qmin"]) # qg[i] = qmin
end

function solve_pf_PVPQ(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_pf_PVPQ; kwargs...)
end

"specification of the formulation agnostic Power Flow model"
function build_pf_PVPQ(pm::AbstractPowerModel)
    variable_bus_voltage(pm, bounded = false)
    variable_gen_power(pm, bounded = false)
    variable_dcline_power(pm, bounded = false)

    for i in ids(pm, :branch)
        expression_branch_power_ohms_yt_from(pm, i)
        expression_branch_power_ohms_yt_to(pm, i)
    end

    constraint_model_voltage(pm)

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        constraint_theta_ref(pm, i)
        constraint_voltage_magnitude_setpoint(pm, i)

        # if multiple generators, fix power generation degeneracies
        if length(ref(pm, :bus_gens, i)) > 1
            for j in collect(ref(pm, :bus_gens, i))[2:end]
                constraint_gen_setpoint_active(pm, j)
                constraint_gen_setpoint_reactive(pm, j)
            end
        end
    end

    ### changed by Lola ###
    # normally this is PV Bus constraints -> constrain pg and vm
    # now, constrain pg and constrain qg of there is a qg violation
    tolerance = 1e-6
    for (i,bus) in ref(pm, :bus)
        constraint_power_balance(pm, i)

        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            for j in ref(pm, :bus_gens, i)
                constraint_gen_setpoint_active(pm, j) # pg[i] = pg
                gen = ref(pm, nw = nw_id_default, :gen, j)

                # if there is a reactive power violation, constrain qg
                if gen["qg"] > gen["qmax"] || isapprox(gen["qg"], gen["qmax"], atol=tolerance) 
                    constraint_gen_setpoint_reactive_max(pm, j) # constraint from top of file
                elseif gen["qmin"] > gen["qg"] || isapprox(gen["qg"], gen["qmin"], atol=tolerance)
                    constraint_gen_setpoint_reactive_min(pm, j) # constraint from top of file
                # else, just make it a pv bus
                else
                    # this assumes inactive generators are filtered out of bus_gens
                    @assert bus["bus_type"] == 2
                    constraint_voltage_magnitude_setpoint(pm, i) # vm[i] = vm
                    # print("voltage enforced")
                end
            end
        end
        
    end
    ### until here ###

    for (i,dcline) in ref(pm, :dcline)
        #constraint_dcline_power_losses(pm, i) not needed, active power flow fully defined by dc line setpoints
        constraint_dcline_setpoint_active(pm, i)

        f_bus = ref(pm, :bus)[dcline["f_bus"]]
        if f_bus["bus_type"] == 1
            constraint_voltage_magnitude_setpoint(pm, f_bus["index"])
        end

        t_bus = ref(pm, :bus)[dcline["t_bus"]]
        if t_bus["bus_type"] == 1
            constraint_voltage_magnitude_setpoint(pm, t_bus["index"])
        end
    end
end

function adjust_PVPQ(network_data_init, nb_iter)
    network_data = deepcopy(network_data_init)
    PF_res1 = solve_pf_PVPQ(network_data, ACPPowerModel, solver)
    update_data!(network_data,PF_res1["solution"]) # update_data! form powermodels.jl (network and ["solution"] have same format)
    global flag_adjust = false 
    for i in 1:nb_iter # check for 8 iterations
        for (g,gen) in network_data["gen"]
            if network_data["bus"]["$(gen["gen_bus"])"]["bus_type"] != 3
                if gen["qmin"] > PF_res1["solution"]["gen"]["$g"]["qg"] || gen["qmax"] < PF_res1["solution"]["gen"]["$g"]["qg"] # check for violation
                    println("PVPQ iterations: ", i)
                    global flag_adjust = true # if there is a violation flag to true
                    
                end
            end
        end
        if flag_adjust # if flag is true, resolve and put flag to false
            PF_res1 = solve_pf_PVPQ(network_data, ACPPowerModel, solver)
            update_data!(network_data,PF_res1["solution"])
            global flag_adjust = false
        else 
            break
        end

    end
    return PF_res1
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