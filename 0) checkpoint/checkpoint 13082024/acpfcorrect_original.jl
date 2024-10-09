
using PowerModels
using JuMP

"
Functions related to creating an ACPF correction layer

"


# we're instantiating a QCRM model: form -> wr

# function for acpf correction 
function build_acpf_correction(pm::AbstractPowerModel)
    variable_bus_voltage(pm) # -> form -> wr
    variable_gen_power(pm) # -> core -> variable
    variable_branch_power(pm) # -> core -> variable
    variable_dcline_power(pm)

    # CHANGE OBJECTIVE
    # objective_min_fuel_and_flow_cost(pm)
    objective_min_pg_cmd_distance(pm)
    #println(objective_min_pg_cmd_distance(pm))

    constraint_model_voltage(pm) # -> form -> wr

    # add gen cmd constraints
    for i in ids(pm, :gen)
        constraint_gen_power_real_command(pm, i)
    end
    

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i) # -> form -> wr
    end

    for i in ids(pm, :bus)
        constraint_power_balance(pm, i) # -> form -> acp
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i) # -> form -> acp
        constraint_ohms_yt_to(pm, i) # -> form -> acp

        constraint_voltage_angle_difference(pm, i) # -> form -> acr

        constraint_thermal_limit_from(pm, i) # -> core -> constraint
        constraint_thermal_limit_to(pm, i) # -> core -> constraint
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i) # -> core -> constraint
    end
    
    #println(pm.model)
end


function objective_min_pg_cmd_distance(pm::AbstractPowerModel; kwargs...)
    pg_ratio = PowerModels.var(pm)[:pg_ratio] = Dict{Int, Any}()

    for (i, gen) in ref(pm, :gen) # Iterate over generators in the network
        pg = PowerModels.var(pm, :pg, i)  # Access the pg variable for generator i
        pg_cmd = gen["pg_cmd"]   # Access the pg_cmd parameter for generator i

        pg_ratio[i] = (1 - pg / pg_cmd)^2
    end


    # Construct the objective expression
    objective_expr = sum(pg_ratio[i] for i in keys(pg_ratio))

    # Set the objective in the JuMP model
    JuMP.@objective(pm.model, Min, objective_expr)

    return objective_expr
end



function constraint_gen_power_real_command(pm::AbstractPowerModel, f_idx)
    pg = PowerModels.var(pm, :pg, f_idx)
    gen_data = ref(pm, :gen)[f_idx]  # Access generator data

    pg_cmd = gen_data["pg_cmd"]
    alpha = gen_data["alpha"]

    if gen_data["pmax"] == 0 # if there is no active power, only add the upper bound on the command constraint
        # JuMP.@constraint(pm.model, (1 - alpha) * pg_cmd <= pg)
        JuMP.@constraint(pm.model, pg <= (1 + alpha) * pg_cmd)
    else
        JuMP.@constraint(pm.model, (1 - alpha) * pg_cmd <= pg)
        JuMP.@constraint(pm.model, pg <= (1 + alpha) * pg_cmd)
    end

end


function push_generator_setpoints(case)
    alpha = 0.15*ones(length(case["gen"]))

    for (gen_key, gen_value) in case["gen"]
        index = parse(Int, gen_key)  # Convert the key to integer index if needed

        # pg_cmd_value = samples[1]["pg_cmd"][index]  # Access corresponding pg_cmd value
        pg_cmd_value = case["gen"]["$(index)"]["pg"]  # add pg_cmd key, copied from Pg
        if pg_cmd_value == 0 # to avoid dividing by zero in objective
            pg_cmd_value = 1e-4
        end
        gen_value["pg_cmd"] = pg_cmd_value  # Assign pg_cmd value to "pg_cmd" key in gen_value

        alpha_value = alpha[index]
        gen_value["alpha"] = alpha_value 
    end


end


# function push_load_setpoints(case, samples)
    
#     for i in 1:length(case["load"])
#         #case["load"]["$(i)"]["pd"] = samples[1]["pd_cmd"][i]
#         case["load"]["$(i)"]["pd"] = samples[1]["pd_cmd"][i]
#     end
# end



"
Functions related to instantiating QC relaxation with variables x = [Pg, Pd], first step, x = [Pg, vj, Pd]

"


# function to instantiate QCRM model 
function build_opf_load(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_load_power(pm) # extra added variable
    variable_branch_power(pm)
    variable_dcline_power(pm)

    # objective_min_fuel_and_flow_cost(pm) # removed objective from build_opf_load

    constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_power_balance_loadvar(pm, i) # modify function with Pd variable in it
    end

    for i in ids(pm, :branch)
        constraint_ohms_yt_from(pm, i)
        constraint_ohms_yt_to(pm, i)

        constraint_voltage_angle_difference(pm, i)

        constraint_thermal_limit_from(pm, i) # -> core -> constraint
        constraint_thermal_limit_to(pm, i) # -> core -> constraint
    end

    for i in ids(pm, :dcline)
        constraint_dcline_power_losses(pm, i) # -> core -> constraint
    end
end


function push_pf(case)
    load = case["load"]
    nb_loads = length(load)

    for (load_key, load_value) in load
        index = parse(Int, load_key)

        pd = load["$(index)"]["pd"]
        qd = load["$(index)"]["qd"]

        sd = sqrt(pd^2 + qd^2)
        pf = pd/sd

        load_value["pf"] = pf
    end

end


"generates variables for both `active` and `reactive` load"
function variable_load_power(pm::AbstractPowerModel; kwargs...)
    variable_load_power_real(pm; kwargs...)
    load_power_imaginary_constant_pf(pm; kwargs...) # add variable for reactive power, but constrain it using pf and variable active power
    # variable_load_power_imaginary(pm; kwargs...) # for now, only add Pd as variable
end


"variable: `pd[j]` for `j` in `load`"
function variable_load_power_real(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pd = PowerModels.var(pm, nw)[:pd] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_pd", 
        start = comp_start_value(ref(pm, nw, :load, i), "pd_start")
    ) # 


    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(pd[i], 0.95*load["pd"]) # lower bound on active power load is 0.6*Pd
            JuMP.set_upper_bound(pd[i], 1.05*load["pd"]) # upper bound is double nominal load
        end
    end

    report && sol_component_value(pm, nw, :load, :pd, ids(pm, nw, :load), pd)
end

function load_power_imaginary_constant_pf(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qd = PowerModels.var(pm, nw)[:qd] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_qd",
        start = comp_start_value(ref(pm, nw, :load, i), "qd_start")
    )

    pd = PowerModels.var(pm, nw)[:pd]

    if bounded
        for (i, load) in ref(pm, nw, :load)
            pf = load["pf"]  # Power factor
            pd_val = pd[i]   # Active power from the first function
            # Apparent power: sd = pd / pf
            # Reactive power: qd = sqrt(sd^2 - pd^2)
            JuMP.@NLconstraint(pm.model, qd[i]^2 == (pd_val / pf)^2 - pd_val^2)
            JuMP.set_lower_bound(qd[i], 0.0)
        end
    end

    report && sol_component_value(pm, nw, :load, :qd, ids(pm, nw, :load), qd)
end

"variable: `qd[j]` for `j` in `load`"
function variable_load_power_imaginary(pm::AbstractPowerModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qd = PowerModels.var(pm, nw)[:qd] = JuMP.@variable(pm.model,
        [i in ids(pm, nw, :load)], base_name="$(nw)_qd",
        start = comp_start_value(ref(pm, nw, :load, i), "qd_start")
    )

    if bounded
        for (i, load) in ref(pm, nw, :load)
            JuMP.set_lower_bound(qd[i], load["qd"])
            JuMP.set_upper_bound(qd[i], load["qd"])
        end
    end

    report && sol_component_value(pm, nw, :load, :qd, ids(pm, nw, :load), qd)
end



# modified from src -> core -> constraint_template
function constraint_power_balance_loadvar(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance_loadvar(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_loads, bus_pd, bus_qd, bus_gs, bus_bs) # modified
end


# modified from src -> form -> acp
function constraint_power_balance_loadvar(pm::QCRMPowerModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_loads, bus_pd, bus_qd, bus_gs, bus_bs) # modified
    vm   = PowerModels.var(pm, n, :vm, i)
    p    = get(PowerModels.var(pm, n),    :p, Dict()); PowerModels._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(PowerModels.var(pm, n),    :q, Dict()); PowerModels._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(PowerModels.var(pm, n),   :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(PowerModels.var(pm, n),   :qg, Dict()); PowerModels._check_var_keys(qg, bus_gens, "reactive power", "generator")
    pd   = get(PowerModels.var(pm, n),   :pd, Dict()); PowerModels._check_var_keys(pd, bus_loads, "active power", "load") # added
    qd   = get(PowerModels.var(pm, n),   :qd, Dict()); PowerModels._check_var_keys(qd, bus_loads, "reactive power", "load") # added
    ps   = get(PowerModels.var(pm, n),   :ps, Dict()); PowerModels._check_var_keys(ps, bus_storage, "active power", "storage")
    qs   = get(PowerModels.var(pm, n),   :qs, Dict()); PowerModels._check_var_keys(qs, bus_storage, "reactive power", "storage")
    psw  = get(PowerModels.var(pm, n),  :psw, Dict()); PowerModels._check_var_keys(psw, bus_arcs_sw, "active power", "switch")
    qsw  = get(PowerModels.var(pm, n),  :qsw, Dict()); PowerModels._check_var_keys(qsw, bus_arcs_sw, "reactive power", "switch")
    p_dc = get(PowerModels.var(pm, n), :p_dc, Dict()); PowerModels._check_var_keys(p_dc, bus_arcs_dc, "active power", "dcline")
    q_dc = get(PowerModels.var(pm, n), :q_dc, Dict()); PowerModels._check_var_keys(q_dc, bus_arcs_dc, "reactive power", "dcline")

    


    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - sum(pd_var[d] for d in bus_loads) # modified
        - sum(pd for (i,pd) in bus_pd) # modified
        - sum(gs for (i,gs) in bus_gs)*vm^2
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - sum(qd_var[d] for d in bus_loads) # modified
        - sum(qd for (i,qd) in bus_qd) # modified
        + sum(bs for (i,bs) in bus_bs)*vm^2
    )

    if PowerModels._IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end
