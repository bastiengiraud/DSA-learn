function build_custom_opf_DTU(data::Dict{String,Any},Pgen_final::Matrix{Float64},Qgen_final::Matrix{Float64},model=Model())


    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    # Add zeros to turn linear objective functions into quadratic ones
    # so that additional parameter checks are not required
    PowerModels.standardize_cost_terms!(data, order=2)

    # use build_ref to filter out inactive components
    ref = PowerModels.build_ref(data)[:it][pm_it_sym][:nw][nw_id_default]

    # ref contains all the relevant system parameters needed to build the OPF model
    # When we introduce constraints and variable bounds below, we use the parameters in ref.

    gen_index=collect(keys(ref[:gen]))
    bus_index=collect(keys(ref[:bus]))

    #NOW WE START DEFINING VARIABLES AND CONSTRAINTS

    # Add voltage angles va for each bus
    @variable(model, va[i in keys(ref[:bus])])
    #[i in keys(ref[:bus])] adds one `va` variable for each bus in the network

    # Add voltage magnitudes vm for each bus
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)
    # note: this vairable also includes the voltage magnitude limits and a starting value

    # Add active power generation variable pg for each generator (including limits)
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])

    # Add reactive power generation variable qg for each generator (including limits)
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    # Add power flow variables p to represent the active power flow for each branch
    @variable(model, -Inf <= p[(l,i,j) in ref[:arcs]] <= Inf)

    # Add power flow variables q to represent the reactive power flow for each branch
    @variable(model, -Inf <= q[(l,i,j) in ref[:arcs]] <= Inf)


    #????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
    # Add power flow limits for each branch. Why we only use rate_a??
    for arc in ref[:arcs]
        (l,i,j) = arc
        branch = ref[:branch][l]
        if haskey(branch, "rate_a")
            JuMP.set_lower_bound(p[arc], -branch["rate_a"])
            JuMP.set_upper_bound(p[arc],  branch["rate_a"])
            JuMP.set_lower_bound(q[arc], -branch["rate_a"])
            JuMP.set_upper_bound(q[arc],  branch["rate_a"])
        end
    end

    # Add power flow variables p_dc to represent the active power flow for each HVDC line
    @variable(model, p_dc[a in ref[:arcs_dc]])
    # Add power flow variables q_dc to represent the reactive power flow at each HVDC terminal
    @variable(model, q_dc[a in ref[:arcs_dc]])

    for (l,dcline) in ref[:dcline]
        f_idx = (l, dcline["f_bus"], dcline["t_bus"])
        t_idx = (l, dcline["t_bus"], dcline["f_bus"])

        JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
        JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])
        JuMP.set_lower_bound(q_dc[f_idx], dcline["qminf"])
        JuMP.set_upper_bound(q_dc[f_idx], dcline["qmaxf"])

        JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
        JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
        JuMP.set_lower_bound(q_dc[t_idx], dcline["qmint"])
        JuMP.set_upper_bound(q_dc[t_idx], dcline["qmaxt"])
    end

    # index representing which side the HVDC line is starting
    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])


    #can be changed to be adapted to what we need
#=
@NLobjective(model, Max,
    sum(log(sum((pg[i]-Pgen_final[i,j])^2 for (i,gen) in ref[:gen])) for j in 1:size(Pgen_final,2))+
    sum(log(sum((qg[i]-Qgen_final[i,j])^2 for (i,gen) in ref[:gen])) for j in 1:size(Qgen_final,2)))
=#

@NLobjective(model, Max,
    sum(log(sum((pg[i]-Pgen_final[i,j])^2 for (i,gen) in ref[:gen])) for j in 1:size(Pgen_final,2))+
    sum(log(sum((qg[i]-Qgen_final[i,j])^2 for (i,gen) in ref[:gen])) for j in 1:size(Qgen_final,2)))


    #@constraint(model, [k = 1:size(Pgen_final,2)], (sum(pg[i]-Pgen_final[i,k])^2  for (i,gen) in ref[:gen]) >= .1)

    # Fix the voltage angle to zero at the reference bus
    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, va[i] == 0)
    end

    # Nodal power balance constraints
    for (i,bus) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )

        # Reactive power balance at node i
        @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )
    end

    # Branch power flow physics and limit constraints
    for (i,branch) in ref[:branch]

        # Build the from and to variable id of the i-th branch, which is a tuple given by (branch id, from bus, to bus)
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        # the from and to sides of a branch are useful to calculate power losses later


        p_fr = p[f_idx]  # p_fr is a reference to the optimization variable p[f_idx]
        q_fr = q[f_idx]  # q_fr is a reference to the optimization variable q[f_idx]
        p_to = p[t_idx]  # p_to is a reference to the optimization variable p[t_idx]
        q_to = q[t_idx]  # q_to is a reference to the optimization variable q[t_idx]
        # adding constraints to p_fr is equivalent to adding constraints to p[f_idx], and so on

        vm_fr = vm[branch["f_bus"]]  # vm_fr is a reference to the optimization variable vm on the from side of the branch
        vm_to = vm[branch["t_bus"]]  # vm_to is a reference to the optimization variable vm on the to side of the branch
        va_fr = va[branch["f_bus"]]  # va_fr is a reference to the optimization variable va on the from side of the branch
        va_to = va[branch["t_bus"]]  # va_fr is a reference to the optimization variable va on the to side of the branch

        # Compute the branch parameters and transformer ratios from the data
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2
        # tap is assumed to be 1.0 on non-transformer branches


        # AC Line Flow Constraints

        # From side of the branch flow
        @NLconstraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        @NLconstraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        # To side of the branch flow
        @NLconstraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        @NLconstraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Phase Angle Difference Limit
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        if haskey(branch, "rate_a")    #why do we need a haskey?????
            @constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
            @constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
        end
    end


    # HVDC line constraints
    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint

        # Build the from variable id of the i-th HVDC line, which is a tuple given by (hvdc line id, from bus, to bus)
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        # Build the to variable id of the i-th HVDC line, which is a tuple given by (hvdc line id, to bus, from bus)
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        # Constraint defining the power flow and losses over the HVDC line
        @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end


    return model, gen_index, bus_index
    
end