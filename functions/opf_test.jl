# activate path and show active packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.status()

using PowerModels
using Ipopt

# include path for power systems data
data_path = raw"C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/cases/static"
file_path_modified = joinpath(data_path, "pglib_opf_case57_ieee_modified.m")
file_path_original = joinpath(data_path, "pglib_opf_case57_ieee.m")

# initialize data
network_data_mod = PowerModels.parse_file(file_path_modified) 
network_data_ori = PowerModels.parse_file(file_path_original) 

nlines_mod = length(network_data_mod["branch"])
nlines_ori = length(network_data_ori["branch"])

nbus_mod = length(network_data_mod["bus"])
nbus_ori = length(network_data_ori["bus"])

ngens = length(network_data_ori["gen"])

# compare results
result_mod = solve_ac_opf(network_data_mod, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
result_ori = solve_ac_opf(network_data_ori, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

# obtain solution
solution_mod = result_mod["solution"]
solution_ori = result_ori["solution"]

# Initialize loss variables
global ploss_mod = 0
global ploss_ori = 0

# Extract line data and compute losses
for i in 1:nlines_mod
    # Calculate real power loss on this line
    branch_mod = abs( abs(solution_mod["branch"]["$i"]["pf"]) - abs(solution_mod["branch"]["$i"]["pt"]) )
    global ploss_mod += branch_mod
end

# Extract line data and compute losses
for i in 1:nlines_ori
    # Calculate real power loss on this line
    branch_ori = abs( abs(solution_ori["branch"]["$i"]["pf"]) - abs(solution_ori["branch"]["$i"]["pt"]) )
    global ploss_ori += branch_ori
end

for i in 1:nlines_ori
    branch_mod = abs(solution_mod["branch"]["$i"]["pf"]) - abs(solution_mod["branch"]["$i"]["pt"])
    branch_ori = abs(solution_ori["branch"]["$i"]["pf"]) - abs(solution_ori["branch"]["$i"]["pt"])

    if abs(branch_mod - branch_ori) > 0.001
        print("Gaat himmel fout bij: ", i)
    end
end

if abs(ploss_mod - ploss_ori) > 0.001
    for i in (1+nlines_ori):nlines_mod
        branch_mod = abs(solution_mod["branch"]["$i"]["pf"]) - abs(solution_mod["branch"]["$i"]["pt"])
        print("loss over trafo $i: ", branch_mod, "\n")
        print("trafo $i is from ", network_data_mod["branch"]["$i"]["f_bus"], " to ", network_data_mod["branch"]["$i"]["t_bus"], "\n")
    end
end

for i in 1:ngens
    pg_mod = solution_mod["gen"]["$i"]["pg"]
    pg_ori = solution_ori["gen"]["$i"]["pg"]
    qg_mod = solution_mod["gen"]["$i"]["qg"]
    qg_ori = solution_ori["gen"]["$i"]["qg"]

    if abs(pg_mod - pg_ori) > 0.01
        print("fy faen completely different Pg for $i: ", "\n")
        print("Pg mod for $i: ", pg_mod, "\n")
        print("Pg ori for $i: ", pg_ori, "\n")
    end

    if abs(qg_mod - qg_ori) > 0.01
        print("fy faen completely different Qg for $i: ", "\n")
        print("Qg mod for $i: ", qg_mod, "\n")
        print("Qg ori for $i: ", qg_ori, "\n")
    end

end

for i in 1:nbus_ori
    vm_mod = solution_mod["bus"]["$i"]["vm"]
    vm_ori = solution_ori["bus"]["$i"]["vm"]
    va_mod = solution_mod["bus"]["$i"]["va"]
    va_ori = solution_ori["bus"]["$i"]["va"]

    if abs(vm_mod - vm_ori) > 0.01
        print("fy faen completely different Vm for $i: ", i)
    end

    if abs(va_mod - va_ori) > 0.01
        print("fy faen completely different Va for $i: ", i)
    end

end



println("Termination status Modified System: ", result_mod["termination_status"])
println("Termination status Original System: ", result_ori["termination_status"])

println("Objective value Modified System: ", result_mod["objective"])
println("Objective value Original System: ", result_ori["objective"])

println("Total Real Power Loss Modified System [MW]: ", ploss_mod*100)
println("Total Real Power Loss Original System [MW]: ", ploss_ori*100)













