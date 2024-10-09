# Data collection algorithm

"""
Optimization based sampling. Afterwards, a MVND is fitted over the sampled data, and the dataset is enhanced.
Note: you sample only feasible space, so you draw a MVND only over feasible samples.

"""

# specify current path
cd("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code")

# activate path and show active packages
using Pkg
Pkg.activate(".")
Pkg.status()

using PowerModels
using JuMP
using Ipopt
using Plots
PowerModels.silence() 

# include other module_methods
include("build_opf_DTU.jl") 
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/method.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/contingency.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/ssa_module.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/dynamics.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/directed_walk.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/acpfcorrect.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/obbt_lu.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/polytope.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/support.jl")
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/functions/write_dfs.jl")

# import initialization module
include("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/init.jl")
using .Initialize

# initialize 
ndim, bus_gen, gen_index, min_lim, max_lim = dim_and_limits_variables(Initialize.network_data, Initialize.variable_loads)
pm, N, vars, header = instantiate_system_QCRM(Initialize.network_data, Initialize.variable_loads)
pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))
op_header = header_full_op(Initialize.network_data, header, pg_numbers, vm_numbers)

# %% Power Flow Exploration
N = Initialize.opt_samples   # number of power flow solutions maximizing the distance between solutions
case = Initialize.network_data

for i in collect(keys(case["gen"]))
   case["gen"][i]["cost"]=[0,0]     #Assign a 0 to all the terms cost to make the objective function 0. If the system had HVDC lines they would need to be changed to 0 too.
end


# Solving the optimal power flow using the standard function of Power Models
OPF_soln=PowerModels.solve_opf(case, ACPPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))


firststatus = println(OPF_soln["termination_status"]) #answers the solution of the model
solvestatus=[]
solvestatusind=string(OPF_soln["termination_status"])
append!(solvestatus, solvestatusind)



llista_gen=[]
llista_bus=[]

for i in collect(keys(case["bus"]))
   append!(llista_bus, convert(AbstractFloat,case["bus"][i]["bus_i"]))
end

for i in collect(keys(case["gen"]))
   if case["gen"][i]["gen_status"]==1  
      append!(llista_gen, convert(AbstractFloat,case["gen"][i]["gen_bus"]))
   end
   if case["gen"][i]["gen_status"]!=1   # The code is not adapted for unit commitment problems. This is a reminder that the mat file has generators not on.
      throw(DomainError("case data error","There are generators with a state variable different than 0"))
   end
end

Pginicial =[]
for ii in collect(keys(case["gen"]))
   if case["gen"][ii]["gen_status"]==1
      append!(Pginicial, convert(AbstractFloat,OPF_soln["solution"]["gen"][ii]["pg"]))
   end
end

Qginicial =[]
for ii in collect(keys(case["gen"]))
   if case["gen"][ii]["gen_status"]==1
      append!(Qginicial, convert(AbstractFloat,OPF_soln["solution"]["gen"][ii]["qg"]))
   end
end

Vminicial =[]
for ii in collect(keys(case["bus"]))
   append!(Vminicial, convert(AbstractFloat,OPF_soln["solution"]["bus"][ii]["vm"]))
end

Vainicial =[]
for ii in collect(keys(case["bus"]))
   append!(Vainicial, convert(AbstractFloat,OPF_soln["solution"]["bus"][ii]["va"]))
end




Vminter1=hcat(Vminicial,llista_bus)
Vminter2=Vminter1[sortperm(Vminter1[:, 2]), :]
Vm=Vminter2[:,1]

Vainter1=hcat(Vainicial,llista_bus)
Vainter2=Vainter1[sortperm(Vainter1[:, 2]), :]
Va=Vainter2[:,1]

Pginter1=hcat(Pginicial,llista_gen)
Pginter2=Pginter1[sortperm(Pginter1[:, 2]), :]
P_gen_def=Pginter2[:,1]

Qginter1=hcat(Qginicial,llista_gen)
Qginter2=Qginter1[sortperm(Qginter1[:, 2]), :]
Q_gen_def=Qginter2[:,1]


function update_result(Pgen_final, Qgen_final, Voltage_mag_final, Voltage_ang_final, N)
    for i=1:N
        model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
        model = build_custom_opf_DTU(case,Pgen_final[:,:],Qgen_final[:,:],model)

        gen_index=model[2]
        bus_index=model[3]

        result = optimize!(model[1])

        if string(termination_status(model[1])) != "LOCALLY_SOLVED"  && string(termination_status(model[1])) != "ALMOST_LOCALLY_SOLVED"
        #if string(termination_status(model[1])) != "LOCALLY_SOLVED"
            print(termination_status(model[1]), "******************", "\n")
            break   # stop the process as soon as the first not optimal solution is found
        end

        Objective_value = objective_value(model[1])

        Pgen=value.(model[1][:pg])
        Pgen_final_before_sortt=Pgen.data
        C=hcat(Pgen_final_before_sortt,gen_index)
        D=C[sortperm(C[:, 2]), :]
        Pgenn=D[:,1]
    
        Qgen=value.(model[1][:qg])
        Qgen_final_before_sortt=Qgen.data
        C1=hcat(Qgen_final_before_sortt,gen_index)
        D1=C1[sortperm(C1[:, 2]), :]
        Qgenn=D1[:,1]
    
        Pgen_final=hcat(Pgen_final,Pgenn)
        Qgen_final=hcat(Qgen_final,Qgenn)
    
    
        Voltage_mag_before_sort=value.(model[1][:vm])
        E=hcat(Voltage_mag_before_sort,bus_index)
        F=E[sortperm(E[:, 2]), :]
        Voltage_mag=F[:,1]
    
        Voltage_ang_before_sort=value.(model[1][:va])
        E1=hcat(Voltage_ang_before_sort,bus_index)
        F1=E1[sortperm(E1[:, 2]), :]
        Voltage_ang=F1[:,1]
    
        Voltage_mag_final=hcat(Voltage_mag_final,Voltage_mag)
        Voltage_ang_final=hcat(Voltage_ang_final,Voltage_ang)
    
        println("Iteration=$i-----------------", "\n")
        println(Objective_value)
        # Check that the solver terminated without an error
        println("The solver termination status is $(termination_status(model[1]))", "\n")
        append!(solvestatus, string(termination_status(model[1])))
        

        Pgen_final = Matrix{Float64}(Pgen_final)
        Qgen_final = Matrix{Float64}(Qgen_final)
        Voltage_mag_final = Matrix{Float64}(Voltage_mag_final)
        Voltage_ang_final = Matrix{Float64}(Voltage_ang_final)


    end
    return Pgen_final,Qgen_final,Voltage_mag_final,Voltage_ang_final
end
  

# %% Call Data Collection function!
P_gen_def, Q_gen_def, Vm_def, Va_def = update_result(P_gen_def,Q_gen_def,Vm,Va,N) #[10x501], [10x501], [39,501], [39,501]


############### now check feasibility of all samples and do sss analysis etc....

# check AC feasibility of all samples
pf_results_feas = []
load_results_feas = []
op_info_feas = []

pf_results_infeas = []
load_results_infeas = []
op_info_infeas = []

global nb_feasible = 0
global nb_infeasible = 0 
global pvpq_feasible = 0
global initial_feasible = 0
global correct_feasible = 0
tollerance = 1e-3

data_mn = Initialize.network_data

# add contingencies
data_mn["area_gens"] = Dict()
contingencies = []
contingencies_gen = []

# Iterate over the branches in network_data
for (i, branch) in data_mn["branch"]
    i = parse(Int, i)
    if (i in Initialize.contingencies_n1)
        push!(contingencies, (idx=i, label="LINE-$(i)", type="branch"))
    end
end

# add contingencies to tightened network data
data_mn["branch_contingencies"] = contingencies
data_mn["gen_contingencies"] = contingencies_gen

for i in 1:N    
    data_opf_verif = deepcopy(data_mn)

    for g in eachindex(pg_numbers)
        data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = P_gen_def[pg_numbers[g],i] 
        data_opf_verif["gen"]["$(pg_numbers[g])"]["qg"] = Q_gen_def[pg_numbers[g],i] 
    end
    for v in eachindex(vm_numbers)
        data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = Vm_def[vm_numbers[v],i] 
        data_opf_verif["bus"]["$(vm_numbers[v])"]["va"] = Va_def[vm_numbers[v],i] 
    end
   #  for d in eachindex(pd_numbers)
   #      data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = sample_ops[i,length(pg_numbers)+length(vm_numbers)+d]  
   #      pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
   #      pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
   #      sd = pd / pf
   #      qd = sqrt(sd^2 - pd^2)
   #      data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
   #  end

    # dictionary placeholder with OP flags. 1 is feasible, 0 is infeasible
    op_flag = Dict(
        "N0" => 1, # flag for base-case feasibility
        "N1" => 1, # flag for N1 feasibility
        "N1_over_volt" => 0.0, # amount over voltage violation
        "N1_under_volt" => 0.0, # amount under voltage violation
        "N1_flow" => 0.0 # amount flow violation
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

    for i in 0:(length(multinetwork["nw"])-1)
        # https://github.com/lanl-ansi/PowerModels.jl/blob/30e3d392fd95f1c5a81abd248ac3acc9462211b8/src/prob/pf.jl#L286-L292
        PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)) # solves non linear pf (unbounded), constraints not enforced
        initial_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(Initialize.network_data, PF_res0["solution"], tollerance)

        if i == 0
            # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
            if initial_feasibility != true    
                op_flag["N0"] = 0
            end
        end

        if i != 0
            # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # PF_res0["termination_status"] == LOCALLY_SOLVED != true && initial_feasibility != true
            if initial_feasibility != true
                op_flag["N1"] = 0
                op_flag["N1_over_volt"] += vm_vio_over
                op_flag["N1_under_volt"] += vm_vio_under
                op_flag["N1_flow"] += sm_vio
            end
        end

    end        

    # check feasibility
    if op_flag["N0"] == 1 && op_flag["N1"] == 1
        global nb_feasible += 1
        global initial_feasible += 1
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
        global nb_infeasible += 1
       

    end
    

end

print("number of feasible samples: ", nb_feasible, "\n")
    
# get list of feasible operating points x
feasible_ops = []
for i in 1:length(pf_results_feas[:])
    op = get_Full_OP(Initialize.network_data, pf_results_feas[i], load_results_feas[i])
    push!(feasible_ops, op)
end

# get list of infeasible operating points
infeasible_ops = []
for i in 1:length(pf_results_infeas[:])
    op = get_Full_OP(Initialize.network_data, pf_results_infeas[i], load_results_infeas[i])
    push!(infeasible_ops, op)
end


# collect feasible and infeasible ops
global_OPs = vcat(feasible_ops, infeasible_ops)
pf_results_total = vcat(pf_results_feas, pf_results_infeas)
nb_ops = length(global_OPs)

op_info_feasible = vcat(op_info_feas)
op_info_infeasible = vcat(op_info_infeas)

feasible_ops = feasible_ops
infeasible_ops = infeasible_ops

# Create a multivariate normal distribution and sample from this distribution
if Initialize.mvnd_sampling == true
    # start sampling from the multivariate normal distribution
    result_mvnd, elapsed_time_mvnd, memory_mvnd, garbage_mvnd = @timed begin
        sample_mvnd(feasible_ops, Initialize.network_data, data_mn, Initialize.mvnd_samples, Initialize.contingencies_n1)
    end

    feasible_ops_mvnd, pf_results_feas_mvnd, op_info_feas_mvnd, infeasible_ops_mvnd, pf_results_infeas_mvnd, op_info_infeas_mvnd, nb_feasible_mvnd, nb_infeasible_mvnd, pvpq_feasible_mvnd, initial_feasible_mvnd, correct_feasible_mvnd, mvnd_sampling_time = result_mvnd

    # write macros
    # df_macros_mvnd(Initialize.mvnd_sampling, elapsed_time_mvnd, mvnd_sampling_time, memory_mvnd, initial_feasible_mvnd, pvpq_feasible_mvnd, correct_feasible_mvnd, Initialize.directory, Initialize.mvnd_macros_filename)

    # current global OPs
    global_OPs = vcat(feasible_ops, feasible_ops_mvnd, infeasible_ops, infeasible_ops_mvnd)
    pf_results_total = vcat(pf_results_feas, pf_results_feas_mvnd, pf_results_infeas, pf_results_infeas_mvnd)
    
    op_info_feasible = vcat(op_info_feas, op_info_feas_mvnd)
    op_info_infeasible = vcat(op_info_infeas, op_info_infeas_mvnd)

    feasible_ops = vcat(feasible_ops, feasible_ops_mvnd)
    infeasible_ops = vcat(infeasible_ops, infeasible_ops_mvnd)
    
end


# Perform small-signal stability analysis on all samples
if Initialize.sss_analysis == true

    # perform small signal stability analysis for both feasible and infeasible samples
    result_sss, elapsed_time_sss, memory_sss, garbage_sss = @timed begin
        sss_evaluation(Initialize.network_data, global_OPs, pg_numbers, vm_numbers, pd_numbers, Initialize.dir_dynamics, Initialize.case_name)
    end

    # obtain the damping ratio, distance to imaginary axis and eigenvalues for all samples
    total_damp, total_dist, total_eigen = result_sss
end

# create dataframes
df_DW_f, df_DW_i = dfs_init(feasible_ops, infeasible_ops, op_info_feasible, op_info_infeasible)

# dataframes including sss analysis
if Initialize.sss_analysis == true
    df_DW_f, df_DW_i = dfs_sss(df_DW_f, df_DW_i, total_damp, total_dist)
end

# Check if column names feasible&infeasible dataframes match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

# Combine the DataFrames vertically
df_DW = vcat(df_DW_f, df_DW_i)

# write dataframe to CSV
output_path_ops = joinpath(Initialize.directory, Initialize.opt_dataset_filename)
CSV.write(output_path_ops, df_DW)  

summary_result_opt(Initialize.mvnd_sampling)


