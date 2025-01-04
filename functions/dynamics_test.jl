# activate path and show active packages
using Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.status()

using PowerModels
using Ipopt
using Sundials
using Plots
using PowerSystems
using PowerSimulationsDynamics
using PowerSystemCaseBuilder

gr()

include("ssa_module.jl")
include("dynamics.jl")
include("method.jl")
include("support.jl")


# include path for power systems data
data_path = joinpath(dirname(@__DIR__), "cases/static/")
file_path = joinpath(data_path, "pglib_opf_case39_epri.m") # pglib_opf_case39_epri
#file_path = "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/1) Smal signal stability/SSA_module/Code_SSAmodule/WSCC_9_bus.raw"
file_path = joinpath(data_path, "pglib_opf_case162_ieee_dtc.m") #  IEEE118_v32.raw, pglib_opf_case162_ieee_dtc, pglib_opf_case14_ieee
#file_path =                  "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/cases/240busWECC_2018_PSS33.raw"

# generator data
dir_dynamics =              joinpath(dirname(@__DIR__), "cases/dynamic/")
case_name =                 "case162"
#file_dyn =                  "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/cases/dynamic/IEEE 39 bus.dyr"
#file_dyn =                  "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/cases/240busWECC_2018_PSS.dyr"

# initialize data
network_data = PowerModels.parse_file(file_path) 
#scale_all_demand(network_data, 0.8)

total_active_demand = sum(load["pd"] for load in values(network_data["load"]))
println("Total Active Demand (P): $total_active_demand MW")

# compare results
result = solve_ac_opf(network_data, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

# obtain solution
solution = result["solution"]
update_data!(network_data, solution)

# perform small signal stability assessment
syst = create_system(network_data) 
machine_data_dict = load_machine_data(dir_dynamics, case_name)
#syst = System(file_path, file_dyn)
construct_dynamic_model(syst, machine_data_dict) # own function from dynamics.jl
#add_dyn_injectors!(syst, file_dyn)
generators = get_components(x -> true, Generator, syst) # , "generator-34-1"
my_thermal_gen = get_component(ThermalStandard, syst, "gen-1") 
#avr = my_thermal_gen.dynamic_injector.avr # get_Ka(avr) etc...


# transform loads to constant impedance
for l in get_components(PowerSystems.StandardLoad, syst)
    transform_load_to_constant_impedance(l)
end

# define the time span
time_span = (0.0, 30.0)

# add a load step
loads = get_components(c -> c isa ElectricLoad, ElectricLoad, syst)
load_names = map(c -> get_name(c), loads)
l_device = get_component(ElectricLoad, syst, "bus12") # PowerSystems.ElectricLoad 14: bus5, 162: bus12
l_change = PowerSimulationsDynamics.LoadChange(1.0, l_device, :P_ref, 0.09) # from 0.5152
l_trip = PowerSimulationsDynamics.LoadTrip(1.0, l_device)

# define simulation
sim_studied = Simulation(ResidualModel, syst , pwd(), time_span)#, l_change) 
x0 = read_initial_conditions(sim_studied)
# setpoints = get_setpoints(sim_studied)
#show_states_initial_value(sim_studied)

# execute simulation and obtain small signal analysis results
res_studied = small_signal_analysis(sim_studied)
execute!(sim_studied, IDA(), dtmax = 0.02, saveat = 0.02, enable_progress_bar = false)
distance, lambda_stability = dist_to_imaxis(res_studied.eigenvalues)
damping, damping_stability = min_damping(res_studied.eigenvalues)
eigenvalues = res_studied.eigenvalues

# plot the results
results = read_results(sim_studied)
p = plot(xlabel = "time", ylabel = "rotor angle [rad]", xlims = (0, 30), legend = :bottomright)

# Loop over all generators and plot their rotor angles
for i in 1:length(network_data["gen"])
    gen_label = "gen-$i"  # Generator label "generator-$(29+i)-1"
    #gen_label = "generator-$(29+i)-1"  # Generator label "generator-$(29+i)-1"
    angle = get_state_series(results, (gen_label, :δ))  # Get rotor angle for each generator
    plot!(angle, label = gen_label, show = true)  # Add to the same plot
end

# Display the plot
display(p)

# Initialize the plot
p2 = plot(xlabel = "time", ylabel = "frequency [Hz]", xlims = (0, 30), legend = :bottomright)

# get base frequency
nominal_frequency_hz = 60  # Or 60 depending on your system
omega_base = 2 * pi * nominal_frequency_hz

# Loop over all generators and plot their rotor angles
for i in 1:length(network_data["gen"])
    gen_label = "gen-$i" # "generator-$(29+i)-1"
    #gen_label = "generator-$(29+i)-1" # "generator-$(29+i)-1"
    
    # Get the state series for omega (assuming :ω is the symbol for omega)
    # The get_state_series likely returns a tuple (time, omega_pu)
    time_series, omega_pu = get_state_series(results, (gen_label, :ω))  # Split the tuple

    # Convert omega from pu to Hz
    frequency_hz = (omega_pu * omega_base) / (2 * pi)

    # Plot the converted frequency in Hz
    plot!(time_series, frequency_hz, label = gen_label, show = true)  # Use time_series as the x-axis
end


# Display the plot
display(p2)
print(damping)


for i in 1:length(network_data["gen"])
    gen_label = "gen-$i"  # Generator label "generator-$(29+i)-1"
    #gen_label = "generator-$(29+i)-1"  # Generator label "generator-$(29+i)-1"
    time_step, angle = get_state_series(results, (gen_label, :δ))  # Get rotor angle for each generator
    if minimum(angle) <= -1
        println("Generator", (network_data["gen"]["$i"]["gen_bus"]) , " has a rotor angle lower than -1")
    end
end

