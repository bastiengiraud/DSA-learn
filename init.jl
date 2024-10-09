###########################################################################
# This module contains all data necessary to initialize the main script.
# These are all the design parameters you can modify for each case.
###########################################################################

module Initialize

using PowerModels
include("functions/acpfcorrect.jl")
include("functions/support.jl")
include("functions/method.jl")

export data_path, file_path, dir_dynamics, case_name, directory, dataset_filename, flows_filename, init_macros_filename, mvnd_macros_filename, network_data, 
mvnd_sampling, contingency_analysis, sss_analysis, directed_walks, variable_loads, hyperplanes, polytope_samples, mvnd_samples, stopping_iteration,
stopping_percentage, contingencies_inf, contingencies_n1, eigenvalues_dws, k_max, k_max_HIC, nominal_load, lhc_dataset_filename, lhc_flows_filename, lhc_samples, opt_dataset_filename,
opt_samples, imp_dataset_filename, imp_flows_filename, lhc_imp_samples, nb_imp_samples, stability_boundary, stability_margin

project_directory =         "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/"
case_number =               "39"

# include path for power systems steady state data
data_path =                 "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/cases/static/"
file_path =                 joinpath(data_path, "pglib_opf_case39_epri.m") # pglib_opf_case39_epri.m pglib_opf_case162_ieee_dtc.m

#file_path = "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/1) Smal signal stability/SSA_module/Code_SSAmodule/WSCC_9_bus.raw"

# file directory for dynamic data
dir_dynamics =              "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/cases/dynamic/"
case_name =                 "case$(case_number)" #

# file directory for storing datasets
directory =                 "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/output/case$(case_number)/datasets/"
dataset_filename =          "ops3.csv"
flows_filename =            "flows3.csv"
init_macros_filename =      "macros3.csv"
mvnd_macros_filename =      "macros_mvnd3.csv"

# LHC sampling initialization
lhc_dataset_filename =      "lhc_ops.csv"
lhc_flows_filename =        "lhc_flows.csv"
lhc_samples =               10000

# optimization based sampling
opt_dataset_filename =      "opt_ops.csv"
opt_samples =               500

# importance sampling
imp_dataset_filename =      "imp_ops2.csv"
imp_flows_filename =        "imp_flows2.csv"
lhc_imp_samples =           2500
nb_imp_samples =            10000

# initialize data
network_data = PowerModels.parse_file(file_path) 
push_load_pf(network_data) # add the powerfactor of the loads to the dataset

# specify sampling methods
sss_analysis =              true
directed_walks =            true
mvnd_sampling =             true
variable_loads =            keys(Dict{Int64, Any}())  # 4 => nothing add the INDICES of the loadvars'

# specify number of hyperplanes and samples in polytope and mvnd 
hyperplanes =               25
polytope_samples =          500
mvnd_samples =              2000

# specify stopping criteria for hyperplane generation
stopping_iteration =        20
stopping_percentage =       0.05

# specify infeasible contingencies, determine which line outages you want to consider for N-1 security
contingencies_n1 =          []
contingencies_inf =         [] 

# define stability boundary 
stability_boundary =        0.03
stability_margin =          0.0025

# specify for how many of the eigenvalues closest to the imaginary axis you want to do DWs
eigenvalues_dws =           25
k_max =                     30
k_max_HIC =                 15

# set load profile
nominal_load = true
scale_all_demand(network_data, 0.8) # scale nominal load profile with constant

if nominal_load != true
    nb_load = length(network_data["load"])
    sample_array = gen_samples(100, nb_load, 0.6, 1)
    feasible, infeasible, non_feasible_index, feasible_index = OPF_feasible_samples(sample_array, network_data)
    sample_array_filtered = [sample_array[i,:] for i in eachindex(sample_array[:,1]) if !(i in non_feasible_index)]
    #demand_profile = sample_array_filtered[2] 
    #update_all_demand_reactive(network_data, demand_profile)
    update_all_demand(network_data, sample_array_filtered[1])
end

end




