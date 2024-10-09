#Base.compilecache(Base.identify_package("GR"))
cd("C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code")

using Pkg
Pkg.activate(".")
Pkg.status()

using CSV
using DataFrames
using Plots
using StatsPlots
using PowerModels
using Ipopt
using JuMP
using JSON
gr()

include("acpfcorrect.jl")
include("obbt.jl")


# include path for power systems data
data_path = raw"C:/Users/bagir/Documents/1) Projects/2) Datasets/pglib-opf"
file_path = joinpath(data_path, "pglib_opf_case30_as.m")

# initialize data
network_data1 = PowerModels.parse_file(file_path) 
push_pf(network_data1) # add the powerfactor of the loads to the dataset

network_data2 = PowerModels.parse_file(file_path) 
push_pf(network_data2) # add the powerfactor of the loads to the dataset

# plotting 
dir_plot_image = "C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/output/case30/images/"
dir_plot_data = "C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/output/case30/datasets/"
dir_plot_obbt = "C:/Users/bagir/Documents/1) Projects/2) Datasets/2) Datasets code/output/case30/obbt/"

# Specified directories
obbt_uload_path = joinpath(dir_plot_obbt, "obbt_uload.json")
obbt_fload_path = joinpath(dir_plot_obbt, "obbt_fload.json")

stats_uncertain_path = joinpath(dir_plot_obbt, "stats_uload.json")
stats_fixed_path = joinpath(dir_plot_obbt, "stats_fload.json")


# # obbt with load uncertainty
# datatight_uload, stats_uload, mod_u = solve_obbt_load_opf!(network_data1, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0), max_iter=3, model_constructor=QCRMPowerModel)

# # standard obbt
# datatight_fload, stats_fload, mod_f = solve_obbt_mod_opf!(network_data2, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0), max_iter=3, model_constructor=QCRMPowerModel)

# # Write the obbt data to JSON
# open(obbt_uload_path, "w") do f
#     JSON.print(f, datatight_uload)
# end

# open(obbt_fload_path, "w") do f
#     JSON.print(f, datatight_fload)
# end

# # Write stats data to JSON
# open(stats_uncertain_path, "w") do f
#     JSON.print(f, stats_uload)
# end

# open(stats_fixed_path, "w") do f
#     JSON.print(f, stats_fload)
# end



# Read the obbt data from JSON
datatight_uload = JSON.parsefile(obbt_uload_path)
datatight_fload = JSON.parsefile(obbt_fload_path)

# Read stats data from JSON
stats_uload = JSON.parsefile(stats_uncertain_path)
stats_fload = JSON.parsefile(stats_fixed_path)


# original vmin and vmax values
vmins_original = [network_data1["bus"][key]["vmin"] for key in keys(network_data1["bus"])]
vmaxs_original = [network_data1["bus"][key]["vmax"] for key in keys(network_data1["bus"])]
vrange_original = vmaxs_original - vmins_original

# uncertain load tightened vmin and vmax values
vmins_uncertain = [datatight_uload["bus"][key]["vmin"] for key in keys(datatight_uload["bus"])]
vmaxs_uncertain = [datatight_uload["bus"][key]["vmax"] for key in keys(datatight_uload["bus"])]
vrange_uncertain = vmaxs_uncertain - vmins_uncertain

# uncertain load tightened vmin and vmax values
vmins_fixed = [datatight_fload["bus"][key]["vmin"] for key in keys(datatight_fload["bus"])]
vmaxs_fixed = [datatight_fload["bus"][key]["vmax"] for key in keys(datatight_fload["bus"])]
vrange_fixed = vmaxs_fixed - vmins_fixed

# volume reductions
vbounds_reduction_uncertain = (vrange_uncertain./vrange_original)*100
vbounds_reduction_fixed = (vrange_fixed./vrange_original)*100



# Prepare data for plotting
df1 = DataFrame(
    value = vcat(vbounds_reduction_uncertain, vbounds_reduction_fixed),
    type = repeat(["Uncertain", "Fixed"], inner=[length(vbounds_reduction_uncertain)])
)

# Create box plots
p1 = @df df1 boxplot(:type, :value, legend=false, xlabel="Condition", ylabel="Reduction (%)", title="Vm reduction")



# original vmin and vmax values
angmins_original = [network_data1["branch"][key]["angmin"] for key in keys(network_data1["branch"])]
angmaxs_original = [network_data1["branch"][key]["angmax"] for key in keys(network_data1["branch"])]
angrange_original = angmaxs_original - angmins_original

# uncertain load tightened vmin and vmax values
angmins_uncertain = [datatight_uload["branch"][key]["angmin"] for key in keys(datatight_uload["branch"])]
angmaxs_uncertain = [datatight_uload["branch"][key]["angmax"] for key in keys(datatight_uload["branch"])]
angrange_uncertain = angmaxs_uncertain - angmins_uncertain

# uncertain load tightened vmin and vmax values
angmins_fixed = [datatight_fload["branch"][key]["angmin"] for key in keys(datatight_fload["branch"])]
angmaxs_fixed = [datatight_fload["branch"][key]["angmax"] for key in keys(datatight_fload["branch"])]
angrange_fixed = angmaxs_fixed - angmins_fixed

# volume reductions
angbounds_reduction_uncertain = (angrange_uncertain./angrange_original)*100
angbounds_reduction_fixed = (angrange_fixed./angrange_original)*100


# Prepare data for plotting
df2 = DataFrame(
    value = vcat(angbounds_reduction_uncertain, angbounds_reduction_fixed),
    type = repeat(["Uncertain", "Fixed"], inner=[length(angbounds_reduction_uncertain)])
)

# Create box plots
p2 = @df df2 boxplot(:type, :value, legend=false, xlabel="Condition", ylabel="Reduction (%)", title="Va reduction")


plot(p1, p2, layout=(1, 2))



image_name = "obbt_volume_reduction.png"
image_path = joinpath(dir_plot_image, image_name)
savefig(image_path)

