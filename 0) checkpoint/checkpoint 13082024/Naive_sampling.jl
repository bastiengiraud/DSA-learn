include("/Users/lolacharles/Downloads/Thesis/Code/Exported/module_methods.jl")
include("/Users/lolacharles/Downloads/Thesis/Code/Dummy/New_module.jl")
include("/Users/lolacharles/Downloads/Thesis/Code/N-1 tests/N_1_methods.jl")

system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case3_lmbd.m"
#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case5_pjm.m"
#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case14_ieee.m"
#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case39_epri.m"
#system_data = "/Users/lolacharles/Downloads/Thesis/Code/TestCases/pglib_opf_case57_ieee.m"

network_data = PowerModels.parse_file(system_data)

ndim, bus_gen, gen_index, min_lim, max_lim = dim_and_limits_variables(network_data)
pm, N, vars, header = instantiate_system_QCRM(network_data)
pg_numbers, vm_numbers = extract_number_and_type(vcat(header[1]))
nb_sample_ops = 1000

function gen_samples_vectors_naive(n_samples, n_dimensions, level_min, level_max)
    scaling_list = [(level_min[i], level_max[i]) for i in eachindex(level_min)]
    plan, _ = LHCoptim(n_samples, n_dimensions, 2)
    scaled_plan = scaleLHC(plan, scaling_list)
    return scaled_plan
end


global sample_ops = gen_samples_vectors_naive(nb_sample_ops,ndim,min_lim,max_lim)
pf_results = []
global nb_feasible =0 
for i in 1:nb_sample_ops      #      nb_samples          #length(optimal_setpoints[:,1]) 
    #println("$i______________________")
    data_opf_verif = deepcopy(network_data)
    
    for g in eachindex(pg_numbers)
        data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = sample_ops[i,g] #optimal_setpoints[i][g] #in_samples[g,i]
    end
    for v in eachindex(vm_numbers)
        data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = sample_ops[i,length(pg_numbers)+v] #optimal_setpoints[i][length(pg_numbers)+v] #in_samples[length(pg_numbers)+v,i]
    end
    
    PF_res1 = adjust_PVPQ(data_opf_verif, 8)
    #PF_res1 = solve_pf(data_opf_verif, ACPPowerModel, solver)
    if PF_res1["termination_status"] == LOCALLY_SOLVED
        global nb_feasible += 1
        push!(pf_results, PF_res1)
    else 
        println(PF_res1["termination_status"] )
    end

end


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

for i in 1:length(pf_results)
    vm_vio_over, vm_vio_under = check_vm_violations(network_data, pf_results[i]["solution"])
    push!(over_array, vm_vio_over)
    push!(under_array, vm_vio_under)
    pg_vio, qg_vio = check_pg_pq_violations(network_data, pf_results[i]["solution"])
    sm_vio = check_flow_violations(network_data, pf_results[i]["solution"])
    push!(pg_array, pg_vio)
    push!(qg_array, qg_vio)
    push!(sm_array, sm_vio)

end
threshold_qg = 10^-06
threshold_sm = 10^-06

nb_pg, index_pg = find_zeros_and_indexes(pg_array,threshold_qg)
nb_qg, index_qg = find_zeros_and_indexes(qg_array,threshold_qg)
nb_sm, index_sm = find_zeros_and_indexes(sm_array,10^-06)
nb_ovi, index_ovi = find_zeros_and_indexes(over_array,10^-06)
nb_uvi, index_uvi = find_zeros_and_indexes(under_array,10^-06)
common_elements = [x for x in index_qg if (x in index_sm  && x in index_ovi && x in index_uvi && x in index_pg )]
rest_of_the_indices = [x for x in eachindex(pf_results) if !(x in common_elements)]

OPs_Naive = []
for i in common_elements
    op = get_Full_OP(network_data, pf_results[i]["solution"])
    push!(OPs_Naive, op)
end

OPs_notFeas_Naive = []
for i in rest_of_the_indices
    op = get_Full_OP(network_data, pf_results[i]["solution"])
    push!(OPs_notFeas_Naive, op)
end

df_DW_f = DataFrame(hcat(OPs_Naive...)', Symbol.(header[1]))
df_DW_f.Feas = ones(nrow(df_DW_f))

df_DW_i = DataFrame(hcat(OPs_notFeas_Naive...)', Symbol.(header[1]))
df_DW_i.Feas = zeros(nrow(df_DW_i))

# Check if column names match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

# Combine the DataFrames vertically
df_DW = vcat(df_DW_i)

CSV.write("/Users/lolacharles/Downloads/Thesis/Code/Final/Datasets/DS_39bus_ACOPF_Naive_bis.csv", df_DW)

# N-1 step
cont_out = ["1"]
#cont_out_14 = ["1"]
#cont_out = ["5","14", "20","27", "32", "33", "34", "37", "39", "41","46" ]

total_sm_array_naive, total_qg_array_naive, total_over_array_naive, total_under_array_naive , total_pg_array_naive , Nb_N_1_naive,result["solution"]["nw"], multinetwork["nw"] = N_1_step(network_data, OPs_Naive, cont_out)


df = CSV.read("/Users/lolacharles/Downloads/Thesis/Code/Final/DS_39bus_ACOPF.csv", DataFrame)
first_eight_columns = df[:, 1:19]
vector_of_arrays = [collect(row) for row in eachrow(first_eight_columns)]
#update_all_demand_reactive(network_data, Demand_39[1])
total_sm_array, total_qg_array, total_over_array, total_under_array , total_pg_array , Nb_N_1, resul, net = N_1_step(network_data,vector_of_arrays, cont_out)




threshold = 1e-6
indices_below_threshold_naive = []
for (i, sub_vector) in enumerate(total_sm_array_naive)
    #if mean(sub_vector) < 0.5

    if (all_below_threshold(total_sm_array_naive[i], threshold) &&
        all_below_threshold(total_under_array_naive[i], threshold) &&
        all_below_threshold(total_over_array_naive[i], threshold) &&
        #all_below_threshold(total_qg_array_naive[i], threshold) &&
        all_below_threshold(total_pg_array_naive[i], threshold))
        push!(indices_below_threshold_naive, i)
    end
end

file_path = "/Users/lolacharles/Downloads/Thesis/Code/Final/Total_damp_5bus.txt"
# Initialize an empty array to store floats
data_array = Float64[]

# Open the file in read mode
open(file_path) do file
    for line in eachline(file)
        # Convert each line to a float and append to data_array
        push!(data_array, parse(Float64, line))
    end
end
