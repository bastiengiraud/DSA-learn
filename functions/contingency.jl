using LatinHypercubeSampling
using JuMP
using Ipopt
using PowerModels
using PowerModelsAnnex
using PowerModelsSecurityConstrained
using PowerModels

include("acpfcorrect.jl")
include("support.jl")

solver = Ipopt.Optimizer


# generator variables but unbounded
function variable_gen_power_modif(pm::AbstractPowerModel; kwargs...)
    variable_gen_power_real_modif(pm; kwargs...)
    variable_gen_power_imaginary_modif(pm; kwargs...)
end


"variable: `pg[j]` for `j` in `gen`"
function variable_gen_power_real_modif(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pg = var(pm, nw)[:pg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_pg",
        start = comp_start_value(ref(pm, nw, :gen, i), "pg_start")
    )
    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(pg[i], gen["pmin"])
            JuMP.set_upper_bound(pg[i], gen["pmax"])
        end
    else # unbounded generator i.e. min = 0, and no max
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(pg[i], 0.0)
            if gen["pmax"] == 0.0
                JuMP.set_upper_bound(pg[i], gen["pmax"])
            end
        end    
    end

    report && sol_component_value(pm, nw, :gen, :pg, ids(pm, nw, :gen), pg)
end

"variable: `qq[j]` for `j` in `gen`"
function variable_gen_power_imaginary_modif(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qg = var(pm, nw)[:qg] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :gen)], base_name="$(nw)_qg",
        start = comp_start_value(ref(pm, nw, :gen, i), "qg_start")
    )

    if bounded
        for (i, gen) in ref(pm, nw, :gen)
            JuMP.set_lower_bound(qg[i], gen["qmin"])
            JuMP.set_upper_bound(qg[i], gen["qmax"])
        end
    end

    report && sol_component_value(pm, nw, :gen, :qg, ids(pm, nw, :gen), qg)
end


function solve_ac_opf_modif(file, optimizer; kwargs...)
    return solve_opf_modif(file, ACPPowerModel, optimizer; kwargs...)
end

function solve_opf_modif(file, model_type::Type, optimizer; kwargs...)
    return solve_model(file, model_type, optimizer, build_opf_modif; kwargs...)
end

# unbounded opf i.e. voltages, generators and flows can exceed their limits
function build_opf_modif(pm::AbstractPowerModel)
    variable_bus_voltage(pm, bounded = false)
    variable_gen_power_modif(pm, bounded = false)
    variable_branch_power(pm , bounded = false)
    variable_dcline_power(pm , bounded = false)
   
    objective_min_fuel_and_flow_cost(pm)

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i)
        constraint_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i)
    end
end

# modified from https://github.com/lanl-ansi/PowerModelsSecurityConstrained.jl/blob/master/src/core/data.jl. Only one modification
function build_c1_scopf_multinetwork_modif(network::Dict{String,<:Any})

    contingencies = length(network["gen_contingencies"]) + length(network["branch_contingencies"])

    #info(_LOGGER, "building scopf multi-network with $(contingencies+1) networks")

    if contingencies > 0
        mn_data = replicate(network, contingencies)
        base_network = mn_data["nw"]["0"] = deepcopy(mn_data["nw"]["1"])

        for (n, network) in mn_data["nw"]
            if n == "0"
                continue
            end

            for (i,bus) in network["bus"]
                if haskey(bus, "evhi")
                    bus["vmax"] = bus["evhi"]
                end
                if haskey(bus, "evlo")
                    bus["vmin"] = bus["evlo"]
                end
            end

            for (i,branch) in network["branch"]
                if haskey(branch, "rate_c")
                    branch["rate_a"] = branch["rate_c"]
                end
            end
        end

        network_id = 1
        for cont in base_network["gen_contingencies"]
            cont_nw = mn_data["nw"]["$(network_id)"]
            cont_nw["name"] = cont.label
            cont_gen = cont_nw["gen"]["$(cont.idx)"]
            cont_gen["gen_status"] = 0

            gen_buses = Set{Int}()
            for (i,gen) in cont_nw["gen"]
                if gen["gen_status"] != 0
                    push!(gen_buses, gen["gen_bus"])
                end
            end
            cont_nw["gen_buses"] = gen_buses

            # commented out. response gens are gens that are allowed corrective control. we only consider preventive SCOPF.
            network["response_gens"] = Set()
            gen_bus = cont_nw["bus"]["$(cont_gen["gen_bus"])"]
            cont_nw["response_gens"] = [] 
            #cont_nw["response_gens"] = cont_nw["area_gens"][gen_bus["area"]] 

            network_id += 1
        end
        for cont in base_network["branch_contingencies"]
            cont_nw = mn_data["nw"]["$(network_id)"]
            cont_nw["name"] = cont.label
            cont_branch = cont_nw["branch"]["$(cont.idx)"]
            cont_branch["br_status"] = 0

            gen_buses = Set{Int}()
            for (i,gen) in cont_nw["gen"]
                if gen["gen_status"] != 0
                    push!(gen_buses, gen["gen_bus"])
                end
            end
            cont_nw["gen_buses"] = gen_buses

            fr_bus = cont_nw["bus"]["$(cont_branch["f_bus"])"]
            to_bus = cont_nw["bus"]["$(cont_branch["t_bus"])"]

            cont_nw["response_gens"] = Set()
            if haskey(cont_nw["area_gens"], fr_bus["area"])
                cont_nw["response_gens"] = cont_nw["area_gens"][fr_bus["area"]]
            end
            if haskey(network["area_gens"], to_bus["area"])
                cont_nw["response_gens"] = union(cont_nw["response_gens"], cont_nw["area_gens"][to_bus["area"]])
            end

            network_id += 1
        end

    else
        mn_data = replicate(network, 1)
        mn_data["nw"]["0"] = mn_data["nw"]["1"]
        delete!(mn_data["nw"], "1")
    end

    return mn_data
end


#########################################
# scopf with load uncertainty 
# modified from https://github.com/lanl-ansi/PowerModelsSecurityConstrained.jl/blob/master/src/prob/scopf.jl

# enables support for v[1], required for objective_variable_pg_cost when pg is an expression
Base.getindex(v::JuMP.GenericAffExpr, i::Int64) = v

function run_c1_scopf_load(file, model_constructor, solver; kwargs...)
    return solve_model(file, model_constructor, solver, build_c1_scopf_load; multinetwork=true, kwargs...)
end


function build_c1_scopf_load(pm::AbstractPowerModel)
    # base-case network id is 0

    variable_bus_voltage(pm, nw=0)
    variable_gen_power(pm, nw=0)
    variable_branch_power(pm, nw=0)
    variable_load_power(pm, variable_loads, nw=0) # extra added variable

    constraint_model_voltage(pm, nw=0)

    for i in ids(pm, :ref_buses, nw=0)
        constraint_theta_ref(pm, i, nw=0)
    end

    for i in ids(pm, :bus, nw=0)
        constraint_power_balance_loadvar(pm, i, nw=0) # modified power balance equations with load var
    end

    for i in ids(pm, :branch, nw=0)
        constraint_ohms_yt_from(pm, i, nw=0)
        constraint_ohms_yt_to(pm, i, nw=0)

        constraint_voltage_angle_difference(pm, i, nw=0)

        constraint_thermal_limit_from(pm, i, nw=0)
        constraint_thermal_limit_to(pm, i, nw=0)
    end


    contigency_ids = [id for id in nw_ids(pm) if id != 0]
    for nw in contigency_ids
        # unbounded variables in contingencies i.e. allow for post-contingency violations.
        variable_bus_voltage(pm, nw=nw) # changed bounded to true. If this is set to false, it raises an error.
        variable_gen_power(pm, nw=nw) # changed bounded to true. 
        variable_branch_power(pm, nw=nw)
        variable_load_power(pm, variable_loads, nw=nw) # extra added variable

        # we consider preventive scopf. No generator response. 
        # variable_c1_response_delta(pm, nw=nw)

        constraint_model_voltage(pm, nw=nw)

        for i in ids(pm, :ref_buses, nw=nw)
            constraint_theta_ref(pm, i, nw=nw)
        end

        gen_buses = ref(pm, :gen_buses, nw=nw)
        for i in ids(pm, :bus, nw=nw)
            constraint_power_balance_loadvar(pm, i, nw=nw) # modified power balance equations with load var

            # if a bus has active generators, fix the voltage magnitude to the base case
            if i in gen_buses
                # constraint_c1_voltage_magnitude_link(pm, i, nw_1=0, nw_2=nw) # doesn't work for QCRM relaxation. All we want is V0 to be equal to Vc
                constraint_c1_voltage_magnitude_link_custom(pm, 0, nw, i) # modified, linking constraint for voltages
            end
        end

        slack_bus = get_slack_idx_mn(pm)
        response_gens = ref(pm, :response_gens, nw=nw)
        for (i,gen) in ref(pm, :gen, nw=nw)
            pg_base = PowerModels.var(pm, :pg, i, nw=0)

            # Skip the slack bus
            if i == slack_bus
                continue 
            end

            # no response gens, all are linked
            constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)

            # # make sure no gens in response_gens
            # response_gens = []

            # # setup the linear response function or fix value to base case
            # if i in response_gens
            #     constraint_c1_gen_power_real_response(pm, i, nw_1=0, nw_2=nw)
            # else
            #     # preventive SCOPF. Set gen power base case to gen power after contingency. 
            #     constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)
            # end
        end


        for i in ids(pm, :branch, nw=nw)
            constraint_ohms_yt_from(pm, i, nw=nw)
            constraint_ohms_yt_to(pm, i, nw=nw)

            constraint_voltage_angle_difference(pm, i, nw=nw)

            constraint_thermal_limit_from(pm, i, nw=nw)
            constraint_thermal_limit_to(pm, i, nw=nw)
        end
    end


end


function constraint_c1_voltage_magnitude_link_custom(pm::QCRMPowerModel, n_1::Int, n_2::Int, i::Int)
    vm_1 = PowerModels.var(pm, n_1, :vm, i)
    vm_2 = PowerModels.var(pm, n_2, :vm, i)

    JuMP.@constraint(pm.model, vm_1 == vm_2)
end

function constraint_c1_voltage_magnitude_link_custom(pm::ACPPowerModel, n_1::Int, n_2::Int, i::Int)
    vm_1 = PowerModels.var(pm, n_1, :vm, i)
    vm_2 = PowerModels.var(pm, n_2, :vm, i)

    JuMP.@constraint(pm.model, vm_1 == vm_2)
end


#######################################

function run_c1_scopf_modif(file, model_constructor, solver; kwargs...)
    return solve_model(file, model_constructor, solver, build_c1_scopf_modif; multinetwork=true, kwargs...)
end




""
# modified from https://github.com/lanl-ansi/PowerModelsSecurityConstrained.jl/blob/master/src/prob/scopf.jl
function build_c1_scopf_modif(pm::AbstractPowerModel)
    # base-case network id is 0

    variable_bus_voltage(pm, bounded = false, nw=0)
    variable_gen_power(pm, bounded = false,nw=0)
    
    variable_branch_power(pm, bounded = false,nw=0)

    constraint_model_voltage(pm, nw=0)

    for i in ids(pm, :ref_buses, nw=0)
        constraint_theta_ref(pm, i, nw=0)
        constraint_voltage_magnitude_setpoint(pm, i, nw=0)

        # if multiple generators, fix power generation degeneracies
        if length(ref(pm,nw=0, :bus_gens, i)) > 1
            for j in collect(ref(pm, nw=0, :bus_gens, i))[2:end]
                constraint_gen_setpoint_active(pm, j, nw=0)
                constraint_gen_setpoint_reactive(pm, j, nw=0)
            end
        end
    end
     
    tolerance = 1e-6
    for i in ids(pm, :bus, nw=0)
        constraint_power_balance(pm, i, nw=0)
         # PV Bus Constraints
         if length(ref(pm,nw=0, :bus_gens, i)) > 0 && !(i in ids(pm, nw=0,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            for j in ref(pm, :bus_gens, i)
                constraint_gen_setpoint_active(pm, j, nw=0)
                
                gen = ref(pm, nw = 0, :gen, j)
                if gen["qg"] > gen["qmax"] || isapprox(gen["qg"], gen["qmax"], atol=tolerance)
                    constraint_gen_setpoint_reactive(pm, 0, gen["index"], gen["qmax"])
                elseif gen["qmin"] > gen["qg"] || isapprox(gen["qg"], gen["qmin"], atol=tolerance)
                    constraint_gen_setpoint_reactive(pm, 0, gen["index"], gen["qmin"])
                else
                    constraint_voltage_magnitude_setpoint(pm, i, nw=0)
                end

            end
        end
    end

    for i in ids(pm, :branch, nw=0)
        constraint_ohms_yt_from(pm, i, nw=0)
        constraint_ohms_yt_to(pm, i, nw=0)

    end


    contigency_ids = [id for id in nw_ids(pm) if id != 0]
    for nw in contigency_ids
        variable_bus_voltage(pm,bounded = false, nw=nw)
        variable_gen_power(pm, bounded = false, nw=nw)
        variable_branch_power(pm, bounded = false, nw=nw)

       


        constraint_model_voltage(pm, nw=nw)

        for i in ids(pm, :ref_buses, nw=nw)
            constraint_theta_ref(pm, i, nw=nw)
            constraint_voltage_magnitude_setpoint(pm, i, nw=nw)

            #if multiple generators, fix power generation degeneracies
            if length(ref(pm,nw=nw, :bus_gens, i)) > 1
               for j in collect(ref(pm,nw=nw, :bus_gens, i))[2:end]
                   constraint_gen_setpoint_active(pm, j, nw=nw)
                   constraint_gen_setpoint_reactive(pm, j, nw=nw)
               end
            end

        end

        gen_buses = ref(pm, :gen_buses, nw=nw)
        for i in ids(pm, :bus, nw=nw)
            constraint_power_balance(pm, i, nw=nw)

            if length(ref(pm,nw=nw, :bus_gens, i)) > 0 && !(i in ids(pm, nw=nw,:ref_buses))
                # this assumes inactive generators are filtered out of bus_gens
    
                #constraint_voltage_magnitude_setpoint(pm, i, nw=nw)
                for j in ref(pm, :bus_gens, i)
                
                    constraint_gen_setpoint_active(pm, j, nw=nw)
                    
                    gen = ref(pm, nw = nw, :gen, j)
                    if gen["qg"] > gen["qmax"] || isapprox(gen["qg"], gen["qmax"], atol=tolerance)
                        constraint_gen_setpoint_reactive(pm, nw, gen["index"], gen["qmax"])
                    elseif gen["qmin"] > gen["qg"] || isapprox(gen["qg"], gen["qmin"], atol=tolerance)
                        constraint_gen_setpoint_reactive(pm, nw, gen["index"], gen["qmin"])
                    else
                        constraint_voltage_magnitude_setpoint(pm, i, nw=nw)
                    end

    
                    
                end
            end
        end

        for i in ids(pm, :branch, nw=nw)
            constraint_ohms_yt_from(pm, i, nw=nw)
            constraint_ohms_yt_to(pm, i, nw=nw)

    
        end
    end


    ##### Setup Objective #####
    
end


#############################





