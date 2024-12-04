
using Logging
using PowerSimulationsDynamics
using UUIDs
using Base.Filesystem: mkpath

include("support.jl")

global_logger(SimpleLogger(stderr, Logging.Error)) # only display errors

#Compute the distance from the imaginary axis 
function dist_to_imaxis(eigenvalues)
    lambda_stability = true

    #Input = eigenvalues of the system -> Vector{Float64}
    real_eigen = real.(eigenvalues)
    imag_eigen = imag.(eigenvalues)
    if maximum(real_eigen) > 10^-8
        println("This state is not small-signal stable.")
        lambda_stability = false
    else -10^-8 < maximum(real_eigen) < 10^-8
        zero_real_indices = findall(real_eigen .== 0)
        zero_real_and_imag_indices = filter(i -> imag_eigen[i] == 0, zero_real_indices)
        #### only print statement when OP is not small-signal stable
        # if isempty(zero_real_and_imag_indices)
        #     println("Warning : 0 but stable")
        # else
        #     println("Warning : there is an eigenvalue at the origin")
        # end
        real_eigen = [eig_val for eig_val in real_eigen if eig_val <= -10^-8]
    end

    dist = maximum(real_eigen)

    return dist, lambda_stability
end

#Compute the minimum damping of the system
function min_damping(eigenvalues; damping_limit = 2)
    damping_stability = true
    damping_threshold = 0.03

    filtered_eigenval = [eig_val for eig_val in eigenvalues if !(-10^-8 <= real.(eig_val) <= 10^-8)]
    real_part_eigen = real.(filtered_eigenval)
    #real_part_eigen = [eig_val for eig_val in real_part_eigen if !(-10^-8 <= eig_val <= 10^-8)]
    im_part_eigen = imag.(filtered_eigenval)
    damping = []
    for i in eachindex(real_part_eigen)
        temp_damp = -real_part_eigen[i]/sqrt(real_part_eigen[i].^2+im_part_eigen[i].^2)
        if abs(temp_damp) <= damping_limit
            push!(damping,temp_damp)
        end
    end

    min_damping = minimum(unique(damping))
    if min_damping < damping_threshold
        damping_stability = false
    end

    return min_damping, damping_stability
end

#Compute the two small-signal indices
function small_signal_module(full_model)
    # Create a unique temporary directory 
    parent_dir = tempdir()
    unique_id = "tempdir_$(UUIDs.uuid4())"
    temp_dir = mktempdir(parent_dir; prefix="$unique_id", cleanup=true)

    ss_stability = false
    sys = full_model
    time_span = (0.0, 30.0)
    sim_studied = Simulation(ResidualModel, sys , temp_dir, time_span) # , initialize_simulation = true)#, initial_conditions = initial_conditions) # origineel
    read_initial_conditions(sim_studied)
    # setpoints = get_setpoints(sim_studied)
    res_studied = small_signal_analysis(sim_studied)
    distance, lambda_stability = dist_to_imaxis(res_studied.eigenvalues)
    damping, damping_stability = min_damping(res_studied.eigenvalues)
    eigenvalues = res_studied.eigenvalues

    if lambda_stability && damping_stability 
        ss_stability = true
    end

    stability = Dict(
        "stability" => Bool(ss_stability),
        "distance" => distance,
        "damping" => damping,
        "eigenvalues" => eigenvalues
    )

    rm(temp_dir; force=true, recursive=true)

    return stability
end


function sss_evaluation(data_tight_HP, global_OPs, pg_numbers, vm_numbers, pd_numbers, dir_dynamics, case_name)

    total_damp = []
    total_dist = []
    total_eigen = []

    data_build = deepcopy(data_tight_HP)
    machine_data_dict = load_machine_data(dir_dynamics, case_name)

    for i in eachindex(global_OPs)   
        # data_build = deepcopy(data_tight_HP)
        for g in eachindex(pg_numbers)
            data_build["gen"]["$(pg_numbers[g])"]["pg"] = global_OPs[i][g] 
        end
        for v in eachindex(vm_numbers)
            data_build["bus"]["$(vm_numbers[v])"]["vm"] = global_OPs[i][length(pg_numbers)+v] 
        end
        for d in eachindex(pd_numbers)
            var_load_index = pd_numbers[d]
            data_build["load"]["$(pd_numbers[d])"]["pd"] = global_OPs[i][length(pg_numbers)+length(vm_numbers)+var_load_index] 
            pf = data_build["load"]["$(pd_numbers[d])"]["pf"]
            pd = data_build["load"]["$(pd_numbers[d])"]["pd"]
            sd = pd / pf
            qd = sqrt(sd^2 - pd^2)
            data_build["load"]["$(pd_numbers[d])"]["qd"] = qd
        end

        print("Small signal analysis number: ", i, "\n")

        syst = create_system(data_build) # builds up data in temp folder
        #add_dyn_system_basic(syst) 
        construct_dynamic_model(syst, machine_data_dict) # builds up data in temp folder. Own function from dynamics.jl
        #add_dyn_injectors!(syst, dir_dynamics) # dyr_content = readlines(Initialize.dir_dynamics)

        stability = small_signal_module(syst)
        push!(total_damp, stability["damping"])
        push!(total_dist, stability["distance"])
        push!(total_eigen, stability["eigenvalues"])

        if i % 500 == 0 # clear temp folder after every 500 iterations
            clean_temp_files()
            GC.gc()
            #clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")
            continue
        end
    
    end
    return total_damp, total_dist, total_eigen
end 






















