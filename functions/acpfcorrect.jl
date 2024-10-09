
"
Functions related to creating an ACPF correction layer

"


using PowerModels
using JuMP


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
Functions related to instantiating QC relaxation with variables x = [Pg, vj, Pd]

"


# function to instantiate QCRM model 
function build_opf_load(pm::AbstractPowerModel)
    variable_bus_voltage(pm)
    variable_gen_power(pm)
    variable_load_power(pm, variable_loads) # extra added variable
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


function push_load_pf(case)
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
function variable_load_power(pm::AbstractPowerModel, variable_loads; kwargs...)
    variable_load_power_real(pm, variable_loads; kwargs...)
    load_power_imaginary_constant_pf(pm, variable_loads; kwargs...) # add variable for reactive power, but constrain it using pf and variable active power
    # variable_load_power_imaginary(pm; kwargs...) # for now, only add Pd as variable
end


"variable: `pd[j]` for `j` in `load`"
function variable_load_power_real(pm::AbstractPowerModel, variable_loads; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    pd = PowerModels.var(pm, nw)[:pd] = JuMP.@variable(pm.model,
        [i in variable_loads], base_name="$(nw)_pd", 
        start = comp_start_value(ref(pm, nw, :load, i), "pd_start")
    ) # [i in ids(pm, nw, :load)]


    if bounded
        #for (i, load) in ref(pm, nw, :load)
        for i in variable_loads
            load = ref(pm, nw, :load)[i]
            JuMP.set_lower_bound(pd[i], 0.95*load["pd"]) # lower bound on active power load is 0.6*Pd
            JuMP.set_upper_bound(pd[i], 1.05*load["pd"]) # upper bound is double nominal load
        end
    end

    #report && sol_component_value(pm, nw, :load, :pd, ids(pm, nw, :load), pd)
    report && sol_component_value(pm, nw, :load, :pd, variable_loads, pd)
end

function load_power_imaginary_constant_pf(pm::AbstractPowerModel, variable_loads; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    qd = PowerModels.var(pm, nw)[:qd] = JuMP.@variable(pm.model,
        [i in variable_loads], base_name="$(nw)_qd",
        start = comp_start_value(ref(pm, nw, :load, i), "qd_start")
    )

    pd = PowerModels.var(pm, nw)[:pd]

    if bounded
        #for (i, load) in ref(pm, nw, :load)
        for i in variable_loads
            load = ref(pm, nw, :load)[i]
            pf = load["pf"]  # Power factor
            pd_val = pd[i]   # Active power from the first function
            # Apparent power: sd = pd / pf
            # Reactive power: qd = sqrt(sd^2 - pd^2)
            JuMP.@constraint(pm.model, qd[i]^2 == (pd_val / pf)^2 - pd_val^2)
            JuMP.set_lower_bound(qd[i], 0.0)
        end
    end

    #report && sol_component_value(pm, nw, :load, :qd, ids(pm, nw, :load), qd)
    report && sol_component_value(pm, nw, :load, :qd, variable_loads, qd)
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



function bus_loads_vars(pm::AbstractPowerModel, nw::Int, i::Int, variable_loads::Set{Int})
    # Retrieve the list of load IDs connected to the bus
    bus_loads = ref(pm, nw, :bus_loads, i)

    # Create a dictionary for load variables
    bus_pd_vars = Dict(k => PowerModels.var(pm, nw)[:pd][k] for k in bus_loads if k in variable_loads)
    bus_qd_vars = Dict(k => PowerModels.var(pm, nw)[:qd][k] for k in bus_loads if k in variable_loads)

    return bus_pd_vars, bus_qd_vars
end

function bus_loads_pars(pm::AbstractPowerModel, nw::Int, i::Int, variable_loads::Set{Int})
    # Retrieve the list of load IDs connected to the bus
    bus_loads = ref(pm, nw, :bus_loads, i)

    # Create dictionaries for load parameters
    bus_pd_pars = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads if !(k in variable_loads))
    bus_qd_pars = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads if !(k in variable_loads))

    return bus_pd_pars, bus_qd_pars
end



# modified from src -> core -> constraint_template
function constraint_power_balance_loadvar(pm::AbstractPowerModel, i::Int; nw::Int=nw_id_default)
    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_arcs_sw = ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    # bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)
    bus_storage = ref(pm, nw, :bus_storage, i)

    # bus_pd = Dict(k => ref(pm, nw, :load, k, "pd") for k in bus_loads)
    # bus_qd = Dict(k => ref(pm, nw, :load, k, "qd") for k in bus_loads)

    # Retrieve load variables and parameters
    variable_loads_set = Set(variable_loads)
    bus_pd_vars, bus_qd_vars = bus_loads_vars(pm, nw, i, variable_loads_set)
    bus_pd_pars, bus_qd_pars = bus_loads_pars(pm, nw, i, variable_loads_set)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance_loadvar(pm, nw, i, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd_vars, bus_qd_vars, bus_pd_pars, bus_qd_pars, bus_gs, bus_bs) # modified
end


# modified from src -> form -> acp
function constraint_power_balance_loadvar(pm::QCRMPowerModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd_vars, bus_qd_vars, bus_pd_pars, bus_qd_pars, bus_gs, bus_bs) # modified
    vm   = PowerModels.var(pm, n, :vm, i)
    p    = get(PowerModels.var(pm, n),    :p, Dict()); PowerModels._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(PowerModels.var(pm, n),    :q, Dict()); PowerModels._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(PowerModels.var(pm, n),   :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(PowerModels.var(pm, n),   :qg, Dict()); PowerModels._check_var_keys(qg, bus_gens, "reactive power", "generator")
    # pd_var   = get(PowerModels.var(pm, n),   :pd, Dict()); PowerModels._check_var_keys(pd_var, bus_loads, "active power", "load") # added
    # qd_var   = get(PowerModels.var(pm, n),   :qd, Dict()); PowerModels._check_var_keys(qd_var, bus_loads, "reactive power", "load") # added

    # Variables for active and reactive power loads
    pd_vars = bus_pd_vars
    qd_vars = bus_qd_vars

    # Parameters for active and reactive power loads
    pd_pars = bus_pd_pars
    qd_pars = bus_qd_pars

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
        #- sum(pd_var[d] for d in bus_loads) # modified
        #- sum(pd for (i,pd) in bus_pd) # modified
        - sum(pd_vars[d] for d in keys(pd_vars)) 
        - sum(pd_pars[d] for d in keys(pd_pars))
        - sum(gs for (i,gs) in bus_gs)*vm^2
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        #- sum(qd_var[d] for d in bus_loads) # modified
        #- sum(qd for (i,qd) in bus_qd) # modified
        - sum(qd_vars[d] for d in keys(qd_vars)) 
        - sum(qd_pars[d] for d in keys(qd_pars))
        + sum(bs for (i,bs) in bus_bs)*vm^2
    )

    if PowerModels._IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end



# modified from src -> form -> acp
function constraint_power_balance_loadvar(pm::ACPPowerModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd_vars, bus_qd_vars, bus_pd_pars, bus_qd_pars, bus_gs, bus_bs) # modified
    vm   = PowerModels.var(pm, n, :vm, i)
    p    = get(PowerModels.var(pm, n),    :p, Dict()); PowerModels._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(PowerModels.var(pm, n),    :q, Dict()); PowerModels._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(PowerModels.var(pm, n),   :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(PowerModels.var(pm, n),   :qg, Dict()); PowerModels._check_var_keys(qg, bus_gens, "reactive power", "generator")
    # pd_var   = get(PowerModels.var(pm, n),   :pd, Dict()); PowerModels._check_var_keys(pd_var, bus_loads, "active power", "load") # added
    # qd_var   = get(PowerModels.var(pm, n),   :qd, Dict()); PowerModels._check_var_keys(qd_var, bus_loads, "reactive power", "load") # added

    # Variables for active and reactive power loads
    pd_vars = bus_pd_vars
    qd_vars = bus_qd_vars

    # Parameters for active and reactive power loads
    pd_pars = bus_pd_pars
    qd_pars = bus_qd_pars

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
        #- sum(pd_var[d] for d in bus_loads) # modified
        #- sum(pd for (i,pd) in bus_pd) # modified
        - sum(pd_vars[d] for d in keys(pd_vars)) 
        - sum(pd_pars[d] for d in keys(pd_pars))
        - sum(gs for (i,gs) in bus_gs)*vm^2
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        #- sum(qd_var[d] for d in bus_loads) # modified
        #- sum(qd for (i,qd) in bus_qd) # modified
        - sum(qd_vars[d] for d in keys(qd_vars)) 
        - sum(qd_pars[d] for d in keys(qd_pars))
        + sum(bs for (i,bs) in bus_bs)*vm^2
    )

    if PowerModels._IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end


# modified from src -> form -> acp
function constraint_power_balance_loadvar(pm::QCRMPowerModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd_vars, bus_qd_vars, bus_pd_pars, bus_qd_pars, bus_gs, bus_bs) # modified
    vm   = PowerModels.var(pm, n, :vm, i)
    p    = get(PowerModels.var(pm, n),    :p, Dict()); PowerModels._check_var_keys(p, bus_arcs, "active power", "branch")
    q    = get(PowerModels.var(pm, n),    :q, Dict()); PowerModels._check_var_keys(q, bus_arcs, "reactive power", "branch")
    pg   = get(PowerModels.var(pm, n),   :pg, Dict()); PowerModels._check_var_keys(pg, bus_gens, "active power", "generator")
    qg   = get(PowerModels.var(pm, n),   :qg, Dict()); PowerModels._check_var_keys(qg, bus_gens, "reactive power", "generator")
    # pd_var   = get(PowerModels.var(pm, n),   :pd, Dict()); PowerModels._check_var_keys(pd_var, bus_loads, "active power", "load") # added
    # qd_var   = get(PowerModels.var(pm, n),   :qd, Dict()); PowerModels._check_var_keys(qd_var, bus_loads, "reactive power", "load") # added

    # Variables for active and reactive power loads
    pd_vars = bus_pd_vars
    qd_vars = bus_qd_vars

    # Parameters for active and reactive power loads
    pd_pars = bus_pd_pars
    qd_pars = bus_qd_pars

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
        #- sum(pd_var[d] for d in bus_loads) # modified
        #- sum(pd for (i,pd) in bus_pd) # modified
        - sum(pd_vars[d] for d in keys(pd_vars)) 
        - sum(pd_pars[d] for d in keys(pd_pars))
        - sum(gs for (i,gs) in bus_gs)*vm^2
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        #- sum(qd_var[d] for d in bus_loads) # modified
        #- sum(qd for (i,qd) in bus_qd) # modified
        - sum(qd_vars[d] for d in keys(qd_vars)) 
        - sum(qd_pars[d] for d in keys(qd_pars))
        + sum(bs for (i,bs) in bus_bs)*vm^2
    )

    if PowerModels._IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end






"""

Functions to solve the power flow and to switch PV buses to PQ buses


"""


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

# modified from https://github.com/lanl-ansi/PowerModels.jl/blob/master/src/prob/pf.jl
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

    ########## modified
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
    PF_res1 = solve_pf_PVPQ(network_data, ACPPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    update_data!(network_data,PF_res1["solution"]) # update_data! from powermodels.jl (network and ["solution"] have same format)
    global flag_adjust = false 
    for i in 1:nb_iter # check for nb_iter iterations
        for (g,gen) in network_data["gen"]
            if network_data["bus"]["$(gen["gen_bus"])"]["bus_type"] != 3
                if gen["qmin"] > PF_res1["solution"]["gen"]["$g"]["qg"] || gen["qmax"] < PF_res1["solution"]["gen"]["$g"]["qg"] # check for reactive power violation
                    println("PVPQ iterations: ", i)
                    global flag_adjust = true # if there is a violation flag to true
                    
                end
            end
        end
        if flag_adjust # if flag is true, resolve and put flag to false
            PF_res1 = solve_pf_PVPQ(network_data, ACPPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
            update_data!(network_data,PF_res1["solution"])
            global flag_adjust = false
        else 
            break
        end

    end
    return PF_res1
end


