cd("C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code")
# install necessary packages
using Pkg
Pkg.activate(".")
Pkg.status()
# Pkg.update("LatinHypercubeSampling")
# Pkg.update("PowerModelsAnnex")
# Pkg.update("PowerModelsSecurityConstrained")
# Pkg.add("MosekTools")
# Pkg.build("StatsPlots")
# Pkg.add("Suppressor")
# ENV["R_HOME"] = "C:\\Program Files\\R\\R-4.3.2"
# Pkg.build("RCall")
# Pkg.add("Polyhedra")
# Pkg.add("LinearAlgebra")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Statistics")
# Pkg.add("Distributions")
# Pkg.build("Plots")
# Pkg.add("Memento")
# Pkg.add("InfrastructureModels")

# include other scripts
#include("module_methods_pgvm.jl")
#include("module_methods_pgpd.jl")
include("module_methods_pgvmpd.jl")
include("New_module.jl")
include("NewN-1.jl")
include("ssa_module.jl")
include("acpfcorrect.jl")
include("obbt.jl")


# include path for power systems data
data_path = raw"C:/Users/bagir/Documents/1) Projects/2) Datasets/pglib-opf"
file_path = joinpath(data_path, "pglib_opf_case30_as.m")

# file directory
directory = "C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/output/case30/datasets/"

# initialize data
network_data = PowerModels.parse_file(file_path) 
push_pf(network_data) # add the powerfactor of the loads to the dataset
variable_loads = keys(Dict{Int64, Any}(1 => nothing)) # add the INDICES of the loadvars

#---------------------- for fixed load profile ----------
method_lola = false

if method_lola == true
    nb_load = length(network_data["load"])
    sample_array = gen_samples(100, nb_load, 0.8, 1)
    feasible, infeasible, non_feasible_index, feasible_index = OPF_feasible_samples(sample_array, network_data)
    sample_array_filtered = [sample_array[i,:] for i in eachindex(sample_array[:,1]) if !(i in non_feasible_index)]
    #demand_profile = sample_array_filtered[2] 
    #update_all_demand_reactive(network_data, demand_profile)
    update_all_demand(network_data, sample_array_filtered[1])
end

#------------------------- hyperplanes and sampling ----------------------
pm, N, vars, header = instantiate_system_QCRM(network_data, variable_loads) # initialize Powermodels.jl model = pm, get number of variables = N, vars = variables, header = variable names
pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1])) # get the number i.e. bus index from the header names
op_header = header_full_op(network_data, header, pg_numbers, vm_numbers)

# Define feasible and infeasible operating points
global total_OPs = []
global total_Infeas_OPs = []

# start OBBT and seperating hyperplanes
Nb_HP = 20 # number of hyperplanes
Nb_insamples = 1000 # number of samples sampled in polyhedron

result, elapsed_time, memory, garbage = @timed begin
    infeasibility_certif_constraints_check_pgvmpd(network_data, Nb_HP, Nb_insamples, variable_loads)
end

OPs_HP, pf_results_HP, data_tight_HP, OPs_notFeas_HP, common_elements_HP, rest_of_the_indices_HP, nb_feasible, pvpq_feasible, initial_feasible, correct_feasible, volumes, HP_time, sampling_time = result

# create dataframe of feasible OPs
df_DW_f = DataFrame(hcat(OPs_HP...)', Symbol.(op_header[1]))
df_DW_f.Feas = ones(nrow(df_DW_f))

# create dataframe of infeasible OPs
if isempty(OPs_notFeas_HP)
    df_DW_i = copy(df_DW_f) 
    empty!(df_DW_i)
else 
    df_DW_i = DataFrame(hcat(OPs_notFeas_HP...)', Symbol.(op_header[1]))
    df_DW_i.Feas = zeros(nrow(df_DW_i))
end


# Check if column names feasible&infeasible dataframes match
if names(df_DW_f) != names(df_DW_i)
    throw(ArgumentError("The DataFrames have different column names"))
end

# Combine the DataFrames vertically
df_DW = vcat(df_DW_f, df_DW_i)
output_path_ops = joinpath(directory, "ops_pgvmpd2.csv")
CSV.write(output_path_ops, df_DW)  

# make dataframe for run specifications
df = DataFrame(Volumes = volumes)
df.elapsed = fill(elapsed_time, nrow(df))
df.HP = fill(HP_time, nrow(df))
df.sampling = fill(sampling_time, nrow(df))
df.memory = fill(memory, nrow(df))

output_path_macros = joinpath(directory, "macros_pgvmpd2.csv")
CSV.write(output_path_macros, df)  



# #--------------------------------- MVND --------------------------------------------
# # start sampling from the multivariate normal distribution
# Nb_MVD = 100
# MND_ops, MND_ops_Infeas, nb_feasible_mvnd = sample_MVND(OPs_HP, network_data, data_tight_HP, Nb_MVD)

# # create dataframe of feasible OPs from MVND
# df_DW_f_MVND = DataFrame(hcat(MND_ops...)', Symbol.(op_header[1]))
# df_DW_f_MVND.Feas = ones(nrow(df_DW_f_MVND))

# # create dataframe of infeasible OPs from MVND
# if isempty(MND_ops_Infeas)
#     df_DW_i_MVND = copy(df_DW_f_MVND) 
#     empty!(df_DW_i_MVND)
# else 
#     df_DW_i_MVND = DataFrame(hcat(MND_ops_Infeas...)', Symbol.(op_header[1]))
#     df_DW_i_MVND.Feas = zeros(nrow(df_DW_i_MVND))
# end

# # Check if column names match
# if names(df_DW_f) != names(df_DW_i)
#     throw(ArgumentError("The DataFrames have different column names"))
# end

# # Combine the DataFrames vertically
# df_DW = vcat(df_DW_f, df_DW_i, df_DW_f_MVND, df_DW_i_MVND)

# CSV.write("C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/output/datasets/ops_mvnd.csv", df_DW)






#----------------------------------------- n-1 step --------------------------------
# # N-1 step
# cont_out = ["1","3"]
# cont_out = ["1", "13", "14"]
# cont_out = ["5","14", "20","27", "32", "33", "34", "37", "39", "41","46" ]
# df = CSV.read("/Users/lolacharles/Downloads/Thesis/Code/Final/DS_39bus_ACOPF.csv", DataFrame)
# first_eight_columns = df[:, 1:19]
# vector_of_arrays = [collect(row) for row in eachrow(first_eight_columns)]
# #update_all_demand_reactive(network_data, Demand_39[1])
# @time total_sm_array, total_qg_array, total_over_array, total_under_array , total_pg_array , Nb_N_1, resul, net = N_1_step(network_data, OPs_HP, cont_out)

# threshold = 1e-6
# indices_below_threshold = []
# for (i, sub_vector) in enumerate(total_sm_array)
#     #if mean(sub_vector) < 0.5

#     if (all_below_threshold(total_sm_array[i], threshold) &&
#         all_below_threshold(total_under_array[i], threshold) &&
#         all_below_threshold(total_over_array[i], threshold) #&&
#         #all_below_threshold(total_qg_array[i], threshold) &&
#         #all_below_threshold(total_pg_array[i], threshold)
#         )
#         push!(indices_below_threshold, i)
#     end
# end

# fli = filter(arr -> all_below_threshold(arr, 0.001), total_over_array)


# sm_index = [mean(arr) for arr in total_sm_array]
# over_index = [mean(arr) for arr in total_over_array]
# under_index = [mean(arr) for arr in total_under_array]

# # Create DataFrames
# df_DW_f = DataFrame(hcat(OPs_HP...)', Symbol.(header[1]))
# df_DW_f.Feas = ones(nrow(df_DW_f))

# df_DW_i = DataFrame(hcat(OPs_notFeas_HP...)', Symbol.(header[1]))
# df_DW_i.Feas = zeros(nrow(df_DW_i))

# # Add N_1 column based on feasibility
# df_DW_f.SM_index = sm_index[1:length(df_DW_f.Feas)]
# df_DW_i.SM_index = -1 * ones(nrow(df_DW_i))

# df_DW_f.Over_index = over_index[1:length(df_DW_f.Feas)]
# df_DW_i.Over_index = -1 * ones(nrow(df_DW_i))

# df_DW_f.Under_index = under_index[1:length(df_DW_f.Feas)]
# df_DW_i.Under_index= -1 * ones(nrow(df_DW_i))


# df_DW_f.N_1 = under_index[1:length(df_DW_f.Feas)]
# df_DW_i.N_1 = -1 * ones(nrow(df_DW_i))

# # Combine the DataFrames
# df_DW = vcat(df_DW_f, df_DW_i)

# df_DW.N_1 = zeros(nrow(df_DW))

# # Set N_1 to 1 for feasible rows
# for idx in indices_below_threshold
#     df_DW.N_1[idx] = 1
# end

# # Check if column names match
# if names(df_DW_f) != names(df_DW_i)
#     throw(ArgumentError("The DataFrames have different column names"))
# end

# CSV.write("/Users/lolacharles/Downloads/Thesis/Code/Final/Datasets/DS_14bus_ACOPF_N_1_Demand1_10k.csv", df_DW)






# #------------------------------ SSS --------------------------------
# # SSS step : 
# global_OPs = vcat(OPs_HP, OPs_notFeas_HP)
# index_res = vcat(common_elements_HP, rest_of_the_indices_HP)
# @time total_damp, total_dist, total_eigen = SSS_eval(data_tight_HP, global_OPs, pg_numbers, vm_numbers, pd_numbers)
# writedlm("C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/output/datasets/Damp_try.txt", total_damp)

# #-------- N1 code ----------
# # # Add N_1 column based on feasibility
# # df_DW_f.SM_index = sm_index[1:length(df_DW_f.Feas)]
# # df_DW_i.SM_index = -1 * ones(nrow(df_DW_i))

# # df_DW_f.Over_index = over_index[1:length(df_DW_f.Feas)]
# # df_DW_i.Over_index = -1 * ones(nrow(df_DW_i))

# # df_DW_f.Under_index = under_index[1:length(df_DW_f.Feas)]
# # df_DW_i.Under_index= -1 * ones(nrow(df_DW_i))

# # df_DW_f.N_1 = under_index[1:length(df_DW_f.Feas)]
# # df_DW_i.N_1 = -1 * ones(nrow(df_DW_i))
# #-----------------------

# # Combine the DataFrames
# df_DW = vcat(df_DW_f, df_DW_i)

# #-------------------------
# # df_DW.N_1 = zeros(nrow(df_DW))

# # # Set N_1 to 1 for feasible rows
# # for idx in indices_below_threshold
# #     df_DW.N_1[idx] = 1
# # end
# #------------------------

# df_DW.SSS = total_damp[1:length(df_DW.Feas)]

# # Check if column names match
# if names(df_DW_f) != names(df_DW_i)
#     throw(ArgumentError("The DataFrames have different column names"))
# end

# CSV.write("C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/output/datasets/ops_sss.csv", df_DW)








# #----------------------- Directed walks ---------------------------
# # only do directed walks with OPs closest to stability boundary
# closet_OPS = closest_to_zero_indices(total_damp,2)
# closet_OPS = closest_to_zero_indices_no_abs(total_damp,2)
# cls_OP = global_OPs[closet_OPS]

# distance = [0.045, 0.02, 0.01]
# alpha = [0.1, 0.01, 0.1, 0.01] # step size
# dw_ops = DW_step(data_tight_HP, index_res, closet_OPS,  pf_results_HP, distance, alpha)



