###########################################################################
# This module contains all data necessary to initialize the main script.
# These are all the design parameters you can modify for each case.
###########################################################################

module Initialize

using PowerModels
include("functions/acpfcorrect.jl")
include("functions/support.jl")
include("functions/method.jl")

export data_path, file_path, dir_dynamics, case_name, directory, dataset_filename, flows_filename, pol_macros_filename, mvnd_macros_filename, method_macros, network_data, 
mvnd_sampling, contingency_analysis, sss_analysis, directed_walks, variable_loads, hyperplanes, polytope_samples, mvnd_samples, stopping_iteration,
stopping_percentage, contingencies_inf, contingencies_n1, dw_computation, k_max, k_max_HIC, distance, alpha, nominal_load, lhc_dataset_filename, lhc_flows_filename, 
lhc_macros, lhc_samples, opt_dataset_filename, opt_samples, imp_dataset_filename, imp_flows_filename, imp_macros, lhc_imp_samples, nb_imp_samples, stability_bound, stability_lb, 
stability_ub, temp_folder

# include path for power systems steady state data
case_number =               "162"                                               # this is the number of buses
tag =			            "SD6"                                               # this is simply a unique tag for the dataset file names
data_path =                 joinpath(@__DIR__, "cases/static/")                 # create a file path for the network data
file_path =                 joinpath(data_path, "pglib_opf_case39_epri.m")      # pglib_opf_case39_epri.m pglib_opf_case162_ieee_dtc.m

# file directory for dynamic data
dir_dynamics =              joinpath(@__DIR__, "cases/dynamic/")                # create a file path for the machine parameter data
case_name =                 "case$(case_number)" #

# file directory for storing datasets
directory =                 joinpath(@__DIR__, "output/case$(case_number)/datasets/") 
dataset_filename =          "$(case_number)bus_$(tag)_method_ops.csv"
flows_filename =            "$(case_number)bus_$(tag)_method_flows.csv"
pol_macros_filename =       "$(case_number)bus_$(tag)_method_macros_polytope.csv"
mvnd_macros_filename =      "$(case_number)bus_$(tag)_method_macros_mvnd.csv"
method_macros =             "$(case_number)bus_$(tag)_macros_method.csv"

# LHC sampling initialization
lhc_dataset_filename =      "$(case_number)bus_lhc_ops.csv"
lhc_flows_filename =        "$(case_number)bus_lhc_flows.csv"
lhc_macros =                "$(case_number)bus_macros_lhc.csv"
lhc_samples =               10000

# optimization based sampling
opt_dataset_filename =      "opt_ops.csv"
opt_samples =               500

# importance sampling
imp_dataset_filename =      "$(case_number)bus_imp_ops.csv"
imp_flows_filename =        "$(case_number)bus_imp_flows.csv"
imp_macros =                "$(case_number)bus_macros_imp.csv"
lhc_imp_samples =           2500
nb_imp_samples =            10000

# initialize data
network_data = PowerModels.parse_file(file_path) 
push_load_pf(network_data)                                                          # add the powerfactor of the loads to the dataset

# specify sampling methods. Indicate if you want to do small-signal stability analysis, directed walks and/or importance sampling (mvnd_sampling)
sss_analysis =              true
directed_walks =            true
mvnd_sampling =             false
variable_loads =            keys(Dict{Int64, Any}())                                # if you want to add uncertain loads, specify them like this (4 => nothing) 

# specify number of hyperplanes and samples in polytope and mvnd 
hyperplanes =               400                                                     # this specifies the number of hyperplanes constructed
polytope_samples =          4000                                                    # this specifies the number of initialization points for the directed walks
mvnd_samples =              1

# specify stopping criteria for hyperplane generation
stopping_iteration =        100                                                     # after this many iterations of no volume reduction, the DW algorithm stops
stopping_percentage =       0.05                                                    # this is the volume percentage the polytope needs to be reduced by for the stopping criterion

# specify infeasible contingencies, determine which line outages you want to consider for N-1 security
contingencies_n1 =          []                                                      # if you want to add line contingencies, specifiy the index of the lines
contingencies_inf =         []                                                      # if you know which lines will result in an infeasible opf, you can exclude them here

# define stability boundary 
stability_bound =           0.03 # this is the minimum damping ratio
#stability_margin =          0.002
stability_lb =              0.0275                                                  # this is the lower bound of the HIC region
stability_ub =              0.0325                                                  # and this is the upper bound

# specify for how many of the eigenvalues closest to the imaginary axis you want to do DWs
dw_computation =            "parallel"                                              # if you can do parallel computing of the DWs, choose: "parallel", otherwise: "series"
k_max =                     40                                                      # maximum number of directed walk steps for each initialization point
k_max_HIC =                 15                                                      # maximum number of directed walk steps when the HIC region is reached
distance =                  [0.015, 0.01, 0.005]
alpha =                     [4, 3, 2, 1] # [2, 1.5, 1, 0.5]

# set temporary directory
parent_temp = "/dev/shm"
custom_temp = "tempdir_$(tag)" # make custom name for temp folder
temp_folder = joinpath(parent_temp, custom_temp) # create path for folder

# set load profile
nominal_load = true
scale_all_demand(network_data, 0.8)                                                 # scale nominal load profile with constant

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




