#Base.compilecache(Base.identify_package("GR"))

cd("C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code")

using Pkg
Pkg.activate(".")
Pkg.status()

using CSV
using DataFrames
using Plots
using StatsPlots
gr()

directory = "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/output/case39/"

ops_path = joinpath(directory, "datasets/ops.csv")
macros_path = joinpath(directory, "datasets/macros.csv")
image_name = "images/ops_spread.png"

OPs = CSV.read(ops_path, DataFrame)  
macros = CSV.read(macros_path, DataFrame)  

# Function to count headers containing a specific substring
count_headers_with_substring(df::DataFrame, substring::String) = 
    sum(occursin(substring, name) for name in names(df))

function filter_feas(df::DataFrame)
    return filter(row -> row[:N0] == 1, df)
end

# count number of variables 
count_pg = count_headers_with_substring(OPs, "PG")
count_vm = count_headers_with_substring(OPs, "VM")
count_pd = count_headers_with_substring(OPs, "PD")

# group dataframes to create new ones
OPs = filter_feas(OPs) # function to filter only feasible setpoints
pg = (OPs[!, 1:count_pg])
vm = (OPs[!, count_pg+1:count_pg+count_vm])
pd = (OPs[!, count_pg+count_vm+1:count_pg+count_vm+count_pd])

# collect data
data_pg = [pg[!, col] for col in names(pg)]
data_vm = [vm[!, col] for col in names(vm)]
data_pd = [pd[!, col] for col in names(pd)]


# Create boxplots with consistent styling
plot_pg = boxplot(
    data_pg,
    title="Generators",
    xlabel="Index",
    ylabel="pg [pu]",
    label="",
    legend=false,
    ylim=(0,10) #12
)

plot_vm = boxplot(
    data_vm,
    title="Voltages",
    xlabel="Index",
    ylabel="vm [pu]",
    label="",
    legend=false,
    ylim =(0.85,1.15)
)

plot_pd = boxplot(
    data_pd,
    title="Loads",
    xlabel="Index",
    ylabel="pd [pu]",
    label="",
    legend=false,
    ylim=(0,4.1) #8
)

# Combine the plots into a single layout
plot(
    plot_pg,
    plot_vm,
    plot_pd,
    layout = (1, 3),
    size = (900, 300),  # Adjust size as needed
    legend = false
)



# image_path = joinpath(directory, image_name)
# savefig(image_path)

