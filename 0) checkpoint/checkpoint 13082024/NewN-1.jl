using LatinHypercubeSampling
# using Plots
using JuMP
using Ipopt
using PowerModels
using PowerModelsAnnex
using PowerModelsSecurityConstrained
using PowerModels
# using MosekTools

solver = Ipopt.Optimizer


function build_c1_scopf_multinetwork_modif(network::Dict{String,<:Any})

    contingencies = length(network["gen_contingencies"]) + length(network["branch_contingencies"])

    #info(_LOGGER, "building scopf multi-network with $(contingencies+1) networks")

    if contingencies > 0
        mn_data = replicate(network, contingencies)
        base_network = mn_data["nw"]["0"] = deepcopy(mn_data["nw"]["1"]) # base case i.e. no contingencies is network 0

        for (n, network) in mn_data["nw"] # solve scopf as multi network
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
                    branch["rate_a"] = branch["rate_c"] # use c rating for contingency cases
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

            #network["response_gens"] = Set()
            #gen_bus = cont_nw["bus"]["$(cont_gen["gen_bus"])"]
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
        mn_data = _PM.replicate(network, 1)
        mn_data["nw"]["0"] = mn_data["nw"]["1"]
        delete!(mn_data["nw"], "1")
    end

    return mn_data
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
        vm_vio_over = vm_vio_over/nb_vm
        vm_vio_under = vm_vio_under/nb_vm
        return vm_vio_over, vm_vio_under 
end


function check_pg_pq_violations(network, solution, tollerance)
    pg_vio = 0.0
    qg_vio = 0.0
    for (i,gen) in network["gen"]
        if gen["gen_status"] != 0
            gen_sol = solution["gen"][i]

            if gen_sol["pg"] < gen["pmin"] - tollerance
                pg_vio += gen["pmin"] - gen_sol["pg"]
            end
            if gen_sol["pg"] > gen["pmax"] + tollerance
                pg_vio += gen_sol["pg"] - gen["pmax"]
            end

            if gen_sol["qg"] < gen["qmin"] - tollerance
                qg_vio += abs(gen["qmin"] - gen_sol["qg"])
            end
            if gen_sol["qg"] > gen["qmax"] + tollerance
                qg_vio += gen_sol["qg"] - gen["qmax"]
            end
        end
    end
    return pg_vio, qg_vio

end

function check_flow_violations(network, solution, tollerance)
    sm_vio = 0.0 #NaN
    if haskey(solution, "branch")
        sm_vio = 0.0
        for (i,branch) in network["branch"]
            if branch["br_status"] != 0
                branch_sol = solution["branch"][i]

                s_fr = abs(branch_sol["pf"])
                s_to = abs(branch_sol["pt"])

                if !isnan(branch_sol["qf"]) && !isnan(branch_sol["qt"])
                    s_fr = sqrt(branch_sol["pf"]^2 + branch_sol["qf"]^2)
                    s_to = sqrt(branch_sol["pt"]^2 + branch_sol["qt"]^2)
                end

                rating = branch["rate_a"]

                if s_fr > rating + tollerance
                    sm_vio += s_fr - rating
                end
                if s_to > rating + tollerance
                    sm_vio += s_to - rating
                end
            end
        end
    end
    return sm_vio
end



function check_ac_feasibility(data_tight_tmp, result, tollerance)
    # check for voltage violations
    vm_vio_over, vm_vio_under = check_vm_violations(data_tight_tmp, result["solution"], tollerance)
    
    # check for generator violations
    pg_vio, qg_vio = check_pg_pq_violations(data_tight_tmp, result["solution"], tollerance)
    
    # check for line violations
    sm_vio = check_flow_violations(data_tight_tmp, result["solution"], tollerance)

    if (vm_vio_over + vm_vio_under + pg_vio + qg_vio + sm_vio) == 0
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





### commented everything out below here for now


# using PowerModelsSecurityConstrained
# include("/Users/lolacharles/Downloads/Thesis/Code/Exported/module_methods.jl")
# include("/Users/lolacharles/Downloads/Thesis/Code/Dummy/New_module.jl")
# include("/Users/lolacharles/Downloads/Thesis/Code/N-1 tests/N_1_methods.jl")
# include("/Users/lolacharles/Downloads/Special course/Small-signal_stability/SSA_module/Code_SSAmodule/SSA_module.jl")
# using LinearAlgebra

# pm_network = PowerModels.parse_file("/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case14_ieee.m")

# #pm_network_bis = deepcopy(pm_network)
# #nb_load = length(pm_network["load"])
# #sample_array = gen_samples(2, nb_load, 0.6, 0.7001)
# #feasible, infeasible, non_feasible_index, feasible_index = OPF_feasible_samples(sample_array, pm_network)
# #sample_array_filtered = [sample_array[i,:] for i in eachindex(sample_array[:,1]) if !(i in non_feasible_index)]
# #demand_profile = sample_array_filtered[2] 
# #pdate_all_demand(pm_network, demand_profile)
# iter = 3
# pm_network, stats_tmp = solve_obbt_opf!(pm_network, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0), max_iter=iter, model_constructor=QCLSPowerModel)
# opf_res = solve_opf(pm_network, QCLSPowerModel, solver)
# update_data!(pm_network,opf_res["solution"])
# set_ac_pf_start_values!(pm_network)
# pm_network["branch"]["3"]["br_status"] = 0
# res_ac = solve_ac_pf(pm_network,solver)
# #pm_network["branch"]["3"]["br_status"] = 0
# #res_ac = solve_ac_pf(pm_network,solver) #adjust_PVPQ(pm_network, 8)   #solve_ac_pf(pm_network,solver)
# #update_data!(pm_network,res_ac["solution"])
# #pg_vio, qg_vio = check_pg_pq_violations(pm_network, res_ac["solution"])
# #sm_vio = check_flow_violations(pm_network, res_ac["solution"])
# #branch_res = calc_branch_flow_ac(pm_network)
# #check_flow_violations(pm_network, opf_res["solution"])

# #update_data!(pm_network,opf_res["solution"])
# #pf_result = solve_ac_pf(pm_network, Ipopt.Optimizer)


# pm_network["area_gens"] = Dict()
# pm_network["gen_contingencies"] = []
# pm_network["branch_contingencies"] = []
# contingencies = []
# contingencies_gen = []
# push!(contingencies, (idx=parse(Int,"7"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"22"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"24"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"36"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"43"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"1"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"2"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"3"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"4"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"6"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"8"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"9"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"10"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"11"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"12"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"13"), label="LINE-7", type="branch"))
# #push!(contingencies, (idx=parse(Int,"14"), label="LINE-7", type="branch"))
# push!(contingencies, (idx=parse(Int,"15"), label="LINE-7", type="branch"))
# for (i, branch) in pm_network["branch"]
#     if i in ["5","14", "20","27", "32", "33", "34", "37", "39", "41","46" ]
#         continue
#     else
#         push!(contingencies, (idx=parse(Int,i), label="LINE-$(i)", type="branch"))
#     end
# #    break
# end

# #for (j, gen) in pm_network["gen"]
# #    push!(contingencies_gen, (idx=parse(Int,j), label="GEN-$(j)", type="gen"))
# #end

# pm_network["branch_contingencies"] = contingencies
# pm_network["gen_contingencies"] = contingencies_gen
# multinetwork = build_c1_scopf_multinetwork_modif(pm_network)
# if multinetwork["per_unit"] == true
#     for (n, network) in multinetwork["nw"]
#         multinetwork["nw"]["$n"]["per_unit"] = true
#     end
# end

# result = run_c1_scopf_modif(multinetwork, ACPPowerModel, Ipopt.Optimizer)
# #update_active_reactive_power_data!(multinetwork["nw"]["0"],result["solution"]["nw"]["0"])
# #result = run_c1_scopf_modif(multinetwork, ACPPowerModel, Ipopt.Optimizer)

# for i in 0:(length(multinetwork["nw"])-1)
#     pg_vio, qg_vio = check_pg_pq_violations(multinetwork["nw"]["$i"], result["solution"]["nw"]["$i"])
#     println("qg_vio : ", qg_vio)
#     sm_vio = check_flow_violations(multinetwork["nw"]["$i"], result["solution"]["nw"]["$i"])
#     println("sm_vio : ", sm_vio)
#     vm_under_tmp, vm_over_tmp = check_vm_violations(multinetwork["nw"]["$i"], result["solution"]["nw"]["$i"])
# end

# n_1_res = []
# term_status = []
# vm_over = []
# vm_under = []
# pg_vio = [] 
# qg_vio = []
# br_status = []
# feasible_cont = 0
# infeasible_cont = 0

# for (n, network) in multinetwork["nw"]
#     if n == "0"
#         continue
#     end
#     n_1_res_temp = solve_ac_pf(multinetwork["nw"]["$n"], Ipopt.Optimizer)

#     #ac_mod = instantiate_model(pm_network_bis, ACPPowerModel, build_pf)
#     #n_1_res_temp = optimize_model!(ac_mod, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

#     #n_1_res_temp = solve_pf(multinetwork["nw"]["$n"] , QCLSPowerModel, Ipopt.Optimizer)
#     #println(n_1_res_temp["termination_status"])
#     push!(term_status, n_1_res_temp["termination_status"])
#     if n_1_res_temp["termination_status"] == LOCALLY_SOLVED
#         global feasible_cont +=1 
#     else
#         global infeasible_cont +=1
#     end
#     push!(n_1_res, n_1_res_temp["solution"])

#     vm_under_tmp, vm_over_tmp = check_vm_violations(network, n_1_res_temp["solution"])
#     pg_vio_tmp, qg_vio_tmp = check_pg_pq_violations(network, n_1_res_temp["solution"])
    
#     #for (j,br) in pm_network_bis["branch"]
#     #    exp_from_tmp = expression_branch_power_ohms_yt_from(ac_mod, parse(Int,j))
#     #    exp_to = expression_branch_power_ohms_yt_to(ac_mod, parse(Int,j))
#     #end

#     push!(br_status, multinetwork["nw"]["$n"]["branch"]["3"]["br_status"])
#     push!(vm_over, vm_over_tmp)
#     push!(vm_under, vm_under_tmp)
#     push!(pg_vio, pg_vio_tmp)
#     push!(qg_vio, qg_vio_tmp)
    
# end
# println(term_status)
# println(br_status)
# println(infeasible_cont)
# println(vm_over)
# println(vm_under)
# println(pg_vio)
# println(qg_vio)
# #println(n_1_res)
# #n_1_res = solve_ac_pf(multinetwork["nw"]["1"], Ipopt.Optimizer)
