include("support.jl")
include("method.jl")
include("contingency.jl")
include("acpfcorrect.jl")

# import initialization module
#init_dir = joinpath(dirname(@__DIR__), "init.jl")
#include(init_dir)
#using .Initialize

#k_max = Initialize.k_max
#k_max_HIC = Initialize.k_max_HIC


function closest_to_zero_indices(arr, N::Int)
    # Create an array of tuples with absolute value and original index
    abs_with_index = [(abs(val), idx) for (idx, val) in enumerate(arr)]
    
    # Sort the array by absolute value
    sorted_abs_with_index = sort(abs_with_index, by = x -> x[1])
    
    # Get the first N indices from the sorted array
    indices = [sorted_abs_with_index[i][2] for i in 1:N]
    
    return indices
end

function closest_to_zero_indices_no_abs(arr, N::Int)
    # Create an array of tuples with absolute value and original index
    abs_with_index = [(val, idx) for (idx, val) in enumerate(arr)]
    
    # Sort the array by absolute value
    sorted_abs_with_index = sort(abs_with_index, by = x -> x[1])
    
    # Get the first N indices from the sorted array
    indices = [sorted_abs_with_index[i][2] for i in 1:N]
    
    return indices
end


function closest_to_boundary_indices(arr, N::Int, xmin::Number, xmax::Number)
    # Filter the array to keep values that are outside the boundary (either less than xmin or greater than xmax)
    outside_with_index = [(val, idx) for (idx, val) in enumerate(arr) if val < xmin || val > xmax]
    
    # Compute the absolute distance to the nearest boundary for each value
    distance_to_boundary = [(val, abs(val < xmin ? xmin - val : val - xmax), idx) for (val, idx) in outside_with_index]
    
    # Sort by the absolute distance to the boundary
    sorted_by_distance = sort(distance_to_boundary, by = x -> x[2])
    
    # Get the first N indices from the sorted array
    indices = [sorted_by_distance[i][3] for i in 1:min(N, length(sorted_by_distance))]
    
    return indices
end



function ops_in_stability_boundary(arr, lower::Float64, upper::Float64)
    # Create an array of tuples with value and original index
    values_with_index = [(val, idx) for (idx, val) in enumerate(arr)]
    
    # Filter the array for values within the specified range
    filtered_values_with_index = filter(x -> lower <= x[1] <= upper, values_with_index)
    
    # Extract the indices from the filtered array
    indices = [x[2] for x in filtered_values_with_index]
    
    return indices
end


function perturbation_forward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)

    damping_forward = []
    dist_forward = []
    for_grad = []

    # perturb all generators in positive direction
    for (g, gen) in data_build["gen"]
        if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0 # check if not slack bus and not synchronous condenser
            # Temporarily store the original value of pg to restore later
            original_pg = data_build["gen"]["$g"]["pg"]

            # Perturb pg value
            updated_value = original_pg + perturbation
            if data_build["gen"]["$g"]["pmax"] < updated_value # check if perturbation doesn't violate generator limits
                data_build["gen"]["$g"]["pg"] = data_build["gen"]["$g"]["pmax"]
                # println(file, "Out of the bounds ")
            else
                data_build["gen"]["$g"]["pg"] = updated_value
            end
            
            # construct system and add dynamic components
            sys_studied = create_system(data_build)
            construct_dynamic_model(sys_studied, dir_dynamics, case_name)

            stability = small_signal_module(sys_studied) # check small signal stability of data with updated generator
            push!(damping_forward, stability["damping"]) 
            push!(dist_forward, stability["distance"])

            # Restore the original pg value to avoid modifying data_build permanently
            data_build["gen"]["$g"]["pg"] = original_pg

        end
    end 

    for i in 1:length(damping_forward[:])
        dOP_k = abs(damping_forward[i] - stability_boundary) # new distance to stability boundary
        dOP = abs(current_damping - stability_boundary) # current distance to stability boundary
        gradient = (dOP_k - dOP)/perturbation # compute forward gradient
        push!(for_grad, gradient)
    end

    return for_grad, damping_forward, dist_forward

end


function perturbation_backward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)

    damping_backward = []
    dist_backward = []
    back_grad = []

    # perturb all active generators in negative direction
    for (g, gen) in data_build["gen"]
        if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0
            # Temporarily store the original value of pg to restore later
            original_pg = data_build["gen"]["$g"]["pg"]

            # Perturb pg value
            updated_value = original_pg - perturbation

            if data_build["gen"]["$g"]["pmin"] > updated_value # check if gen limits are not violated
                data_build["gen"]["$g"]["pg"] = data_build["gen"]["$g"]["pmin"]
                # println(file, "Out of the bounds ")
            else
                data_build["gen"]["$g"]["pg"] = updated_value
            end
            
            sys_studied = create_system(data_build)
            construct_dynamic_model(sys_studied, dir_dynamics, case_name)

            stability = small_signal_module(sys_studied) # perform small signal stability analysis with negatively perturbed generator
            push!(damping_backward, stability["damping"]) 
            push!(dist_backward, stability["distance"])

            # Restore the original pg value to avoid modifying data_build permanently
            data_build["gen"]["$g"]["pg"] = original_pg

        end
    end 

    for i in 1:length(damping_backward[:])
        dOP_k = abs(damping_backward[i] - stability_boundary) # new distance to stability boundary
        dOP = abs(current_damping - stability_boundary) # current distance to stability boundary
        gradient = (dOP - dOP_k)/perturbation # compute backward gradient
        push!(back_grad, gradient)
    end

    return back_grad, damping_backward, dist_backward

end


function max_gradients(data_build, for_grad, back_grad)

    max_grad = []
    #for i in 1:length(for_grad[:])
    counter = 1
    for (g, gen) in data_build["gen"]
        if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0

            if abs(for_grad[counter]) > abs(back_grad[counter])
                push!(max_grad, for_grad[counter])
            else
                push!(max_grad, back_grad[counter])
            end
            counter += 1
        end
        
    end

    return max_grad

end


function mask_gens_at_limits(data_build, max_grad)
    " If the generator is at its limit, and the gradient is still directed towards that limit
    set the gradient to zero. Otherwise, the directed walk won't do anything."

    mask = falses(length(max_grad))
    counter = 1
    for (g, gen) in data_build["gen"]
        if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"][g]["pmax"] > 0.0

            if data_build["gen"]["$g"]["pg"] == data_build["gen"]["$g"]["pmax"] && max_grad[counter] < 0
                mask[counter] = true
            elseif data_build["gen"]["$g"]["pg"] == data_build["gen"]["$g"]["pmin"] && max_grad[counter] > 0
                mask[counter] = true
            end
            counter += 1
        end
        
    end

    return mask

end


function dw_step_gen(data_build, epsilon, max_grad)

    ind = 1
    for (g, gen) in data_build["gen"]
        if data_build["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_build["gen"]["$g"]["pmax"] > 0.0
            step_size = epsilon*data_build["gen"]["$g"]["pmax"] # scale step size with generator limits
            data_build["gen"]["$g"]["pg"] -= step_size * max_grad[ind]

            if data_build["gen"]["$g"]["pmax"] < data_build["gen"]["$g"]["pg"] # check if gen limits are not violated
                data_build["gen"]["$g"]["pg"] = data_build["gen"]["$g"]["pmax"]
            end

            if data_build["gen"]["$g"]["pmin"] > data_build["gen"]["$g"]["pg"] # check if gen limits are not violated
                data_build["gen"]["$g"]["pg"] = data_build["gen"]["$g"]["pmin"]
            end

            ind += 1
        end
    end

end

function array_exists(arr_list, target_array, tol=1e-2) # tollerance is 0.01 pu i.e. 1 MW
    for arr in arr_list
        if length(arr) == length(target_array) && all(isapprox(arr[i], target_array[i], atol=tol) for i in 1:length(arr))
            return true
        end
    end
    return false
end

function remove_duplicate_ops!(feasible_ops, infeasible_ops, dw_ops, dw_stability, tol=1e-2)
    i = 1
    while i <= length(dw_ops)
        if array_exists(feasible_ops, dw_ops[i], tol) || array_exists(infeasible_ops, dw_ops[i], tol)
            # Remove both the operating point and its corresponding stability value
            splice!(dw_ops, i)
            splice!(dw_stability, i)
        else
            i += 1
        end
    end
end


# Function to compute the Euclidean distance between two arrays
function euclidean_distance(x::Vector, y::Vector)
    return norm(x - y)
end

# Function to remove arrays that lie within a given radius R from each other
function remove_nearby_arrays(arrays::Vector{Any}, damp_pol_feas::Vector{Any}, R::Float64)
    # Step 1: Create a mask where damp_pol_feas >= 0 (i.e., stable points)
    valid_mask = damp_pol_feas .>= 0
    
    # Step 2: Create a boolean mask to track which arrays should be kept, initialized to valid_mask
    keep = valid_mask

    # Step 3: Iterate through each pair of arrays where damp_pol_feas >= 0
    for i in 1:length(arrays)
        if keep[i]
            for j in (i+1):length(arrays)
                if keep[j] && euclidean_distance(arrays[i], arrays[j]) < R
                    keep[j] = false  # Mark the j-th array to be removed
                end
            end
        end
    end

    # Step 4: Find the indices where the mask is true (arrays to keep)
    indices = findall(keep)

    # Step 5: Filter and return the arrays that are kept, along with their original indices
    return arrays[indices], indices
end



function DW_step(data_tight, feasible_ops_polytope, infeasible_ops_polytope, closest_ops, cls_op, variable_loads, pf_results_prev, distance, alpha, stability_boundary, stability_lower_bound, stability_upper_bound, dir_dynamics, case_name, k_max, k_max_HIC)

    dw_dir = joinpath(@__DIR__, "DW_file_try.txt")
    file = open(dw_dir, "w")

    pm, N, vars, header = instantiate_system_QCRM(data_tight, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))

    perturbation = 1e-6 # perturbation size of generator

    plot_damping = []
    plot_dist = []

    directed_walk_ops = []
    directed_walk_stability = []

    gen_keys = collect(keys(data_tight["gen"]))
    gen_slack = get_slack_idx(data_tight)

    # remove slack bus, remove synchronous condensers
    filtered_gen_keys = filter(key -> key != gen_slack && data_tight["gen"][key]["pmax"] > 0.0, gen_keys)

    # Get the list of active generators
    active_gens = [data_tight["gen"][key] for key in filtered_gen_keys]

    # Compute the number of active generators
    lgt_active_gen = length(active_gens)

    # loop over the selected operating points close to the stability boundary
    index = 1

    # avoid modifying original data
    data_build = deepcopy(data_tight)

    for i in closest_ops 
        # clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

        # update steady state data with setpoint data
        update_data!(data_build, pf_results_prev[i])
        variable_loads = collect(variable_loads)
        for load in variable_loads
            data_build["load"]["$load"]["pd"] = cls_op[index][length(pg_numbers)+length(vm_numbers)+load]
            pf = data_build["load"]["$load"]["pf"]
            pd = data_build["load"]["$load"]["pd"]
            sd = pd / pf
            qd = sqrt(sd^2 - pd^2)
            data_build["load"]["$load"]["qd"] = qd
        end
        index += 1

        # create system, add dynamic elements, perform SSA
        sys_studied = create_system(data_build) # this builds up temp folder
        construct_dynamic_model(sys_studied, dir_dynamics, case_name) # this also builds up temp folder

        stability = small_signal_module(sys_studied)

        println(file, "Initial damping :", stability["damping"])
        println(file, "Initial distance :", stability["distance"])   
        println(file, "Initial eigenvalue :", stability["eigenvalues"])   

        current_damping = stability["damping"]
        damping_array_global = []
        dist_array_global = []

        HIC_reached = false # no DW steps in HIC area

        # maximum on number of directed walk steps
        for DW in 1:k_max

            print("OP number: ", (index - 1),"/", length(closest_ops), ", Directed walk number: $DW __________________", "\n")

            # if finished with DWs in HIC, skip to next OP
            if HIC_reached == true
                break
            end

            # compute damping ratio with perturbed generator setpoints     
            for_grad, damping_forward, dist_forward = perturbation_forward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)
            back_grad, damping_backward, dist_backward = perturbation_backward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)           

            ############# compute step size #######
            dOP_k = abs(current_damping - stability_boundary)
            
            # define step size of directed walks based on distance
            if dOP_k > distance[1]
                epsilon = alpha[1]
            elseif distance[2] < dOP_k  < distance[1]
                epsilon = alpha[2]
            elseif distance[3] < dOP_k < distance[2]
                epsilon = alpha[3] 
            elseif dOP_k < distance[3]
                epsilon = alpha[4] 
            end

            ############ check to which side the step should be taken
            println(file, "Directed walk number: $DW __________________", "\n")
            max_grad = max_gradients(data_build, for_grad, back_grad)
            sP = get_SetPoint(data_build) # get the generator setpoints
            dw_step_gen(data_build, epsilon, max_grad) # take a directed walk step along all dimensions
            new_SP = get_SetPoint(data_build) # get new generator setpoints

            # some print statements to check directed walks
            println(file, "these are the max gradients: ", join(max_grad), ", ") # write the current generator setpoints to file
            #println(file, "these are the forward gradients: ", join(for_grad), ", ") # write the current generator setpoints to file
            #println(file, "these are the backward gradients: ", join(back_grad), ", ") # write the current generator setpoints to file
            println(file, "this is the epsilon: ", join(epsilon), ", ") # write the current generator setpoints to file
            println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
            println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
        
            # check small signal stability of perturbered generator and obtain stability indices
            sys_studied = create_system(data_build)
            construct_dynamic_model(sys_studied, dir_dynamics, case_name)

            stability = small_signal_module(sys_studied)

            # some print statements to check directed walks
            println(file, "damping of new setpoint: ", stability["damping"])
            println(file, "distance of new setpoint: ", stability["distance"])

            push!(plot_damping, stability["damping"])
            push!(plot_dist, stability["distance"])

            OP_tmp = get_Full_OP(data_tight, data_build, data_build) # use data_tight, then solution, then data_build
            
            # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
            if !array_exists(directed_walk_ops, OP_tmp) && !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp)
                push!(directed_walk_ops, OP_tmp)
                push!(directed_walk_stability, (stability["damping"], stability["distance"]))
            else
                print("This OP already exists, not adding it to the dataset.", "\n")
            end
            
            # push!(directed_walk_ops, OP_tmp)
            # push!(directed_walk_stability, (stability["damping"], stability["distance"]))

            current_damping = stability["damping"]

            print("current damping is :", current_damping, "\n")

            if current_damping < 0 # stop directed walks if unstable point
                break
            end

            if (stability_lower_bound < current_damping) && 
                (current_damping < stability_upper_bound)

                ########### sample points around first HIC point
                print("OP number: ", (index - 1),"/", length(closest_ops), ", Sampling around first HIC point __________________", "\n")
                data_surround = deepcopy(data_build)
                surrounding_ops = []

                # take a 1MW step in both directions for every generator
                for (g, gen) in data_surround["gen"]
                    if data_surround["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_surround["gen"]["$g"]["pmax"] > 0.0
                        sP = get_SetPoint(data_surround) # get the generator setpoints
                        data_surround["gen"]["$g"]["pg"] += 0.01 # 1MW step
                        new_SP = get_SetPoint(data_surround) # get new generator setpoints

                        if data_surround["gen"]["$g"]["pmax"] < data_surround["gen"]["$g"]["pg"] # check if gen limits are not violated
                            data_surround["gen"]["$g"]["pg"] = data_surround["gen"]["$g"]["pmax"]
                        end

                        # check small signal stability of perturbered generator and obtain stability indices
                        sys_studied = create_system(data_surround)
                        construct_dynamic_model(sys_studied, dir_dynamics, case_name)
                        stability = small_signal_module(sys_studied)

                        println(file, "Surrounding HIC setpoint up step __________________", "\n")
                        println(file, "damping of surround setpoint: ", stability["damping"])
                        println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
                        println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
                        
                        # get setpoint
                        OP_tmp = get_Full_OP(data_tight, data_surround, data_surround) # use data_tight, then solution, then data_build
            
                        # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
                        if !array_exists(directed_walk_ops, OP_tmp) && !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp)
                            push!(directed_walk_ops, OP_tmp)
                            push!(surrounding_ops, OP_tmp)
                            push!(directed_walk_stability, (stability["damping"], stability["distance"]))
                            print("New HIC OP found, adding it to the dataset.", "\n")
                        else
                            print("This OP already exists, not adding it to the dataset.", "\n")
                        end

                        sP = get_SetPoint(data_surround) # get the generator setpoints
                        data_surround["gen"]["$g"]["pg"] -= 0.02 # 1MW step to original point, and 1MW step backwards
                        new_SP = get_SetPoint(data_surround) # get new generator setpoints

                        if data_surround["gen"]["$g"]["pmin"] > data_surround["gen"]["$g"]["pg"] # check if gen limits are not violated
                            data_surround["gen"]["$g"]["pg"] = data_surround["gen"]["$g"]["pmin"]
                        end

                        # check small signal stability of perturbered generator and obtain stability indices
                        sys_studied = create_system(data_surround)
                        construct_dynamic_model(sys_studied, dir_dynamics, case_name)
                        stability = small_signal_module(sys_studied)

                        println(file, "Surrounding HIC setpoint down step__________________", "\n")
                        println(file, "damping of surround setpoint: ", stability["damping"])
                        println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
                        println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
                
                        # get setpoint
                        OP_tmp = get_Full_OP(data_tight, data_surround, data_surround) # use data_tight, then solution, then data_build
            
                        # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
                        if !array_exists(directed_walk_ops, OP_tmp) && !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp)
                            push!(directed_walk_ops, OP_tmp)
                            push!(surrounding_ops, OP_tmp)
                            push!(directed_walk_stability, (stability["damping"], stability["distance"]))
                            print("New HIC OP found, adding it to the dataset.", "\n")
                        else
                            print("This OP already exists, not adding it to the dataset.", "\n")
                        end

                        # reset to initial value for to evaluate next gen
                        data_surround["gen"]["$g"]["pg"] += 0.01 # 1MW step

                    end
                end


                ############ continue directed walks along a single dimension
                # get direction of steepest gradient descent along a single dimension
                for DW_HIC in 1:k_max_HIC

                    print("OP number: ", (index - 1),"/", length(closest_ops), ", HIC directed walk number: $DW_HIC __________________", "\n")

                    # compute damping ratio with perturbed generator setpoints     
                    for_grad, damping_forward, dist_forward = perturbation_forward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)
                    back_grad, damping_backward, dist_backward = perturbation_backward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)           
        
                    ############# take the smallest step size #######
                    epsilon = alpha[4] 
                    
                    ############ check to which side the step should be taken
                    println(file, "In HIC DW number: $DW_HIC __________________", "\n")
                    max_grad = max_gradients(data_build, for_grad, back_grad)
                    mask = mask_gens_at_limits(data_build, max_grad)
                    max_grad[mask] .= 0 #  don't consider generators at limit

                    idx_max = findmax(abs, max_grad)[2] # get index of largest absolute gradient
                    val_max = max_grad[idx_max] # get value of largest absolute gradient
                    max_grad .= 0
                    max_grad[idx_max] = val_max # keep only largest gradient for DW step

                    sP = get_SetPoint(data_build) # get the generator setpoints
                    dw_step_gen(data_build, epsilon, max_grad) # take a directed walk step along all dimensions
                    new_SP = get_SetPoint(data_build) # get new generator setpoints

                    # some print statements to check directed walks
                    println(file, "these are the gradients: ", join(max_grad), ", ") # write the current generator setpoints to file
                    println(file, "this is the epsilon: ", join(epsilon), ", ") # write the current generator setpoints to file
                    println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
                    println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
                
                    # check small signal stability of perturbered generator and obtain stability indices
                    sys_studied = create_system(data_build)
                    construct_dynamic_model(sys_studied, dir_dynamics, case_name)

                    stability = small_signal_module(sys_studied)

                    # some print statements to check directed walks
                    println(file, "damping of new setpoint: ", stability["damping"])
                    println(file, "distance of new setpoint: ", stability["distance"])
        
                    push!(plot_damping, stability["damping"])
                    push!(plot_dist, stability["distance"])
        
                    OP_tmp = get_Full_OP(data_build, data_build, data_build)

                    # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
                    # do add if it's part of directed_walk_ops but also surrounding_ops
                    if !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp) && 
                        (!array_exists(directed_walk_ops, OP_tmp) || array_exists(surrounding_ops, OP_tmp))
                        push!(directed_walk_ops, OP_tmp)
                        push!(directed_walk_stability, (stability["damping"], stability["distance"]))
                    else
                        print("This OP already exists, not adding it to the dataset.")
                        break # stop if you found existing OP
                    end
        
                    # push!(directed_walk_ops, OP_tmp)
                    # push!(directed_walk_stability, (stability["damping"], stability["distance"]))
        
                    current_damping = stability["damping"]

                    print("current damping is :", current_damping, "\n")

                    if (stability_lower_bound > current_damping) || 
                        (current_damping > stability_upper_bound) # stop if you left the HIC
                        break
                    end
                
                end

                HIC_reached = true

            end

        end
        
    end
    close(file)

    return directed_walk_ops, directed_walk_stability

end



function DW_step_single_op(data_tight, feasible_ops_polytope, infeasible_ops_polytope, op_number, closest_op, cls_op, variable_loads, pf_results_prev, distance, alpha, stability_boundary, stability_lower_bound, stability_upper_bound, dir_dynamics, case_name, k_max, k_max_HIC)

    #dw_dir = joinpath(@__DIR__, "DW_file_try_$(op_number).txt")
    #file = open(dw_dir, "w")

    pm, N, vars, header = instantiate_system_QCRM(data_tight, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))

    perturbation = 1e-6 # perturbation size of generator

    plot_damping = []
    plot_dist = []

    directed_walk_ops = []
    directed_walk_stability = []

    gen_keys = collect(keys(data_tight["gen"]))
    gen_slack = get_slack_idx(data_tight)

    # remove slack bus, remove synchronous condensers
    filtered_gen_keys = filter(key -> key != gen_slack && data_tight["gen"][key]["pmax"] > 0.0, gen_keys)

    # Get the list of active generators
    active_gens = [data_tight["gen"][key] for key in filtered_gen_keys]

    # Compute the number of active generators
    lgt_active_gen = length(active_gens)

    # avoid modifying original data
    data_build = deepcopy(data_tight)

    # clear_temp_folder("C:/Users/bagir/AppData/Local/Temp")

    # update steady state data with setpoint data
    update_data!(data_build, pf_results_prev[closest_op])
    variable_loads = collect(variable_loads)
    for load in variable_loads
        data_build["load"]["$load"]["pd"] = cls_op[length(pg_numbers)+length(vm_numbers)+load]
        pf = data_build["load"]["$load"]["pf"]
        pd = data_build["load"]["$load"]["pd"]
        sd = pd / pf
        qd = sqrt(sd^2 - pd^2)
        data_build["load"]["$load"]["qd"] = qd
    end

    # create system, add dynamic elements, perform SSA
    sys_studied = create_system(data_build) # this builds up temp folder
    construct_dynamic_model(sys_studied, dir_dynamics, case_name) # this also builds up temp folder

    stability = small_signal_module(sys_studied)

    #println(file, "Initial damping :", stability["damping"])
    #println(file, "Initial distance :", stability["distance"])   
    #println(file, "Initial eigenvalue :", stability["eigenvalues"])   

    current_damping = stability["damping"]
    damping_array_global = []
    dist_array_global = []

    HIC_reached = false # no DW steps in HIC area

    # maximum on number of directed walk steps
    for DW in 1:k_max

        print("OP number: ", op_number, ", Directed walk number: $DW ________", current_damping, "\n")

        # if finished with DWs in HIC, skip to next OP
        if HIC_reached == true
            break
        end

        # compute damping ratio with perturbed generator setpoints     
        for_grad, damping_forward, dist_forward = perturbation_forward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)
        back_grad, damping_backward, dist_backward = perturbation_backward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)           

        ############# compute step size #######
        dOP_k = abs(current_damping - stability_boundary)
        
        # define step size of directed walks based on distance
        if dOP_k > distance[1]
            epsilon = alpha[1]
        elseif distance[2] < dOP_k  < distance[1]
            epsilon = alpha[2]
        elseif distance[3] < dOP_k < distance[2]
            epsilon = alpha[3] 
        elseif dOP_k < distance[3]
            epsilon = alpha[4] 
        end

        ############ check to which side the step should be taken
        #println(file, "Directed walk number: $DW __________________", "\n")
        max_grad = max_gradients(data_build, for_grad, back_grad)
        sP = get_SetPoint(data_build) # get the generator setpoints
        dw_step_gen(data_build, epsilon, max_grad) # take a directed walk step along all dimensions
        new_SP = get_SetPoint(data_build) # get new generator setpoints

        # some print statements to check directed walks
        #println(file, "these are the max gradients: ", join(max_grad), ", ") # write the current generator setpoints to file
        #println(file, "these are the forward gradients: ", join(for_grad), ", ") # write the current generator setpoints to file
        #println(file, "these are the backward gradients: ", join(back_grad), ", ") # write the current generator setpoints to file
        #println(file, "this is the epsilon: ", join(epsilon), ", ") # write the current generator setpoints to file
        #println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
        #println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
    
        # check small signal stability of perturbered generator and obtain stability indices
        sys_studied = create_system(data_build)
        construct_dynamic_model(sys_studied, dir_dynamics, case_name)

        stability = small_signal_module(sys_studied)

        # some print statements to check directed walks
        #println(file, "damping of new setpoint: ", stability["damping"])
        #println(file, "distance of new setpoint: ", stability["distance"])

        push!(plot_damping, stability["damping"])
        push!(plot_dist, stability["distance"])

        OP_tmp = get_Full_OP(data_tight, data_build, data_build) # use data_tight, then solution, then data_build
        
        # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
        if !array_exists(directed_walk_ops, OP_tmp) && !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp)
            #push!(directed_walk_ops, OP_tmp)
            #push!(directed_walk_stability, (stability["damping"], stability["distance"]))
            print("Not adding OPs outside of HIC region.", "\n")
        else
            print("This OP already exists, not adding it to the dataset.", "\n")
        end
        
        current_damping = stability["damping"]

        print("current damping is :", current_damping, "\n")
	
	if (0.02 < current_damping) && (current_damping < 0.04)
            #push!(directed_walk_ops, OP_tmp)
	    #push!(directed_walk_stability, (stability["damping"], stability["distance"]))
	    #print("Adding medium HIC point to dataset, current damping: ", current_damping, "\n")
        end

        if current_damping < 0 # stop directed walks if unstable point
            break
        end

        if (stability_lower_bound < current_damping) && 
            (current_damping < stability_upper_bound)

            # add the first point in the HIC region
            push!(directed_walk_ops, OP_tmp)
            push!(directed_walk_stability, (stability["damping"], stability["distance"]))

            ########### sample points around first HIC point
            print("OP number: ", op_number, ", Sampling around first HIC point _________", current_damping, "\n")
            data_surround = deepcopy(data_build)
            surrounding_ops = []

            # take a 1MW step in both directions for every generator
            for (g, gen) in data_surround["gen"]
                if data_surround["bus"]["$(gen["gen_bus"])"]["bus_type"] !=3 && data_surround["gen"]["$g"]["pmax"] > 0.0
                    sP = get_SetPoint(data_surround) # get the generator setpoints
                    data_surround["gen"]["$g"]["pg"] += 0.01 # 1MW step
                    new_SP = get_SetPoint(data_surround) # get new generator setpoints

                    if data_surround["gen"]["$g"]["pmax"] < data_surround["gen"]["$g"]["pg"] # check if gen limits are not violated
                        data_surround["gen"]["$g"]["pg"] = data_surround["gen"]["$g"]["pmax"]
                    end

                    # check small signal stability of perturbered generator and obtain stability indices
                    sys_studied = create_system(data_surround)
                    construct_dynamic_model(sys_studied, dir_dynamics, case_name)
                    stability = small_signal_module(sys_studied)

                    #println(file, "Surrounding HIC setpoint up step __________________", "\n")
                    #println(file, "damping of surround setpoint: ", stability["damping"])
                    #println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
                    #println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
                    
                    # get setpoint
                    OP_tmp = get_Full_OP(data_tight, data_surround, data_surround) # use data_tight, then solution, then data_build
        
                    # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
                    if !array_exists(directed_walk_ops, OP_tmp) && !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp)
                        push!(directed_walk_ops, OP_tmp)
                        push!(surrounding_ops, OP_tmp)
                        push!(directed_walk_stability, (stability["damping"], stability["distance"]))
                        print("New HIC OP found, adding it to the dataset.", "\n")
                    else
                        print("This OP already exists, not adding it to the dataset.", "\n")
                    end

                    sP = get_SetPoint(data_surround) # get the generator setpoints
                    data_surround["gen"]["$g"]["pg"] -= 0.02 # 1MW step to original point, and 1MW step backwards
                    new_SP = get_SetPoint(data_surround) # get new generator setpoints

                    if data_surround["gen"]["$g"]["pmin"] > data_surround["gen"]["$g"]["pg"] # check if gen limits are not violated
                        data_surround["gen"]["$g"]["pg"] = data_surround["gen"]["$g"]["pmin"]
                    end

                    # check small signal stability of perturbered generator and obtain stability indices
                    sys_studied = create_system(data_surround)
                    construct_dynamic_model(sys_studied, dir_dynamics, case_name)
                    stability = small_signal_module(sys_studied)

                    #println(file, "Surrounding HIC setpoint down step__________________", "\n")
                    #println(file, "damping of surround setpoint: ", stability["damping"])
                    #println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
                    #println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
            
                    # get setpoint
                    OP_tmp = get_Full_OP(data_tight, data_surround, data_surround) # use data_tight, then solution, then data_build
        
                    # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
                    if !array_exists(directed_walk_ops, OP_tmp) && !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp)
                        push!(directed_walk_ops, OP_tmp)
                        push!(surrounding_ops, OP_tmp)
                        push!(directed_walk_stability, (stability["damping"], stability["distance"]))
                        print("New HIC OP found, adding it to the dataset.", "\n")
                    else
                        print("This OP already exists, not adding it to the dataset.", "\n")
                    end

                    # reset to initial value for to evaluate next gen
                    data_surround["gen"]["$g"]["pg"] += 0.01 # 1MW step

                end
            end


            ############ continue directed walks along a single dimension
            # get direction of steepest gradient descent along a single dimension
            for DW_HIC in 1:k_max_HIC

                print("OP number: ", op_number, ", HIC directed walk number: $DW_HIC ________", current_damping, "\n")

                # compute damping ratio with perturbed generator setpoints     
                for_grad, damping_forward, dist_forward = perturbation_forward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)
                back_grad, damping_backward, dist_backward = perturbation_backward(data_build, dir_dynamics, case_name, current_damping, perturbation, stability_boundary)           
    
                ############# take the smallest step size #######
                epsilon = alpha[4] 
                
                ############ check to which side the step should be taken
                #println(file, "In HIC DW number: $DW_HIC __________________", "\n")
                max_grad = max_gradients(data_build, for_grad, back_grad)
                mask = mask_gens_at_limits(data_build, max_grad)
                max_grad[mask] .= 0 #  don't consider generators at limit

                idx_max = findmax(abs, max_grad)[2] # get index of largest absolute gradient
                val_max = max_grad[idx_max] # get value of largest absolute gradient
                max_grad .= 0
                max_grad[idx_max] = val_max # keep only largest gradient for DW step

                sP = get_SetPoint(data_build) # get the generator setpoints
                dw_step_gen(data_build, epsilon, max_grad) # take a directed walk step along all dimensions
                new_SP = get_SetPoint(data_build) # get new generator setpoints

                # some print statements to check directed walks
                #println(file, "these are the gradients: ", join(max_grad), ", ") # write the current generator setpoints to file
                #println(file, "this is the epsilon: ", join(epsilon), ", ") # write the current generator setpoints to file
                #println(file, "current generator setpoints: ", join(collect(values(sP)), ", ")) # write the current generator setpoints to file
                #println(file, "updated generator setpoints: ", join(collect(values(new_SP)), ", ")) # write the updated generator setpoints to file
            
                # check small signal stability of perturbered generator and obtain stability indices
                sys_studied = create_system(data_build)
                construct_dynamic_model(sys_studied, dir_dynamics, case_name)

                stability = small_signal_module(sys_studied)

                # some print statements to check directed walks
                #println(file, "damping of new setpoint: ", stability["damping"])
                #println(file, "distance of new setpoint: ", stability["distance"])
    
                push!(plot_damping, stability["damping"])
                push!(plot_dist, stability["distance"])
    
                OP_tmp = get_Full_OP(data_build, data_build, data_build)

                # Push OP_tmp only if it doesn't exist in the list within the specified tolerance
                # do add if it's part of directed_walk_ops but also surrounding_ops
                if !array_exists(feasible_ops_polytope, OP_tmp) && !array_exists(infeasible_ops_polytope, OP_tmp) && 
                    (!array_exists(directed_walk_ops, OP_tmp) || array_exists(surrounding_ops, OP_tmp))
                    push!(directed_walk_ops, OP_tmp)
                    push!(directed_walk_stability, (stability["damping"], stability["distance"]))
                else
                    print("This OP already exists, not adding it to the dataset.")
                    break # stop if you found existing OP
                end
    
                # push!(directed_walk_ops, OP_tmp)
                # push!(directed_walk_stability, (stability["damping"], stability["distance"]))
    
                current_damping = stability["damping"]

                print("current damping is :", current_damping, "\n")

                if (stability_lower_bound > current_damping) || 
                    (current_damping > stability_upper_bound) # stop if you left the HIC
                    break
                end
            
            end

            HIC_reached = true

        end

    end
        
    #close(file)
    GC.gc()

    return directed_walk_ops, directed_walk_stability

end



function dw_ops_feasibility(network_basic, data_tight_tmp, variable_loads, directed_walk_ops, directed_walk_stability, dir_dynamics, case_name)

    pm, N, vars, header = instantiate_system_QCRM(data_tight_tmp, variable_loads)
    pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))

    pf_results_feas = []
    load_results_feas = []
    op_info_feas = []
    op_info_feas_stability = []

    pf_results_infeas = []
    load_results_infeas = []
    op_info_infeas = []
    op_info_infeas_stability = []
    
    nb_feasible_dws = 0
    nb_infeasible_dws = 0
    initial_feasible_dws = 0
    tollerance = 1e-4

    # avoid modifying original data
    data_opf_verif = deepcopy(data_tight_tmp)

    for counter in 1:(length(directed_walk_ops[:])) 
        
        for g in eachindex(pg_numbers)
            data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = directed_walk_ops[counter][g]
        end
        for v in eachindex(vm_numbers)
            data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = directed_walk_ops[counter][length(pg_numbers)+v] 
        end
        for d in eachindex(pd_numbers)
            var_load_index = pd_numbers[d]
            data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = directed_walk_ops[counter][length(pg_numbers)+length(vm_numbers)+var_load_index] 
            pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
            pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
            sd = pd / pf
            qd = sqrt(sd^2 - pd^2)
            data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
        end

        print("Current DW sample number: ", Int(counter), "\n")

        # dictionary placeholder with OP flags. 1 is feasible, 0 is infeasible
        op_flag = Dict(
            "N0" => 1, # flag for base-case feasibility
            "N0P" => 0.0, # active power violation
            "N0Q" => 0.0, # active power violation
            "N0OV" => 0.0, # amount over voltage violation
            "N0UV" => 0.0, # amount under voltage violation
            "N0L" => 0.0, # amount line flow violation
            "N1" => 1, # flag for N1 feasibility
            "N1OV" => 0.0, # amount over voltage violation
            "N1UV" => 0.0, # amount under voltage violation
            "N1L" => 0.0 # amount line flow violation
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
        pg_vio = 0.0
        qg_vio = 0.0

        for i in 0:(length(multinetwork["nw"])-1)
            PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
            update_data!(multinetwork["nw"]["$i"], PF_res0["solution"]) # update data with PF results
            flows0 = calc_branch_flow_ac(multinetwork["nw"]["$i"]) # compute branch flows
            update_data!(multinetwork["nw"]["$i"], flows0) # add branch flows
            update_data!(PF_res0["solution"], flows0) # add branch flows to solution
            initial_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res0["solution"], tollerance)

            print(pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio)
            if i == 0
                # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # initial_feasibility != true
                if initial_feasibility != true 
                    op_flag["N0"] = 0
                    op_flag["N0P"] += pg_vio
                    op_flag["N0Q"] += qg_vio
                    op_flag["N0OV"] += vm_vio_over
                    op_flag["N0UV"] += vm_vio_under
                    op_flag["N0L"] += sm_vio
                end
            end

            if i != 0
                # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # initial_feasibility != true
                if initial_feasibility != true 
                    op_flag["N1"] = 0
                    op_flag["N1OV"] += vm_vio_over
                    op_flag["N1UV"] += vm_vio_under
                    op_flag["N1L"] += sm_vio
                end
            end

        end                

        # check feasibility
        if op_flag["N0"] == 1 && op_flag["N1"] == 1
            nb_feasible_dws += 1
            initial_feasible_dws += 1
            push!(pf_results_feas, PF_res0["solution"])
            push!(load_results_feas, data_opf_verif)
            push!(op_info_feas, op_flag) 
            push!(op_info_feas_stability, directed_walk_stability[counter])
            println("initial status:", PF_res0["termination_status"] , "\n")
            print("initial feasibility: ", initial_feasibility, "\n")
        else 

            vars_new = []

            pm = instantiate_model(multinetwork, ACPPowerModel, build_c1_scopf_load)

            # get generator variables and add to vars
            slack_gen_idx = get_slack_idx_mn(pm)
            for i=1:length(pm.var[:it][:pm][:nw][0][:pg].data)
                if i != slack_gen_idx && data_tight_tmp["gen"]["$i"]["pmax"] > 0.0
                    push!(vars_new, JuMP.variable_by_name(pm.model, string("0_pg[",i,"]")))
                end
            end

            # get voltage magnitude variables and add to vars
            gen_indexes = unique(map(x -> x["gen_bus"], values(pm.data["nw"]["0"]["gen"])))
            for g in gen_indexes
                push!(vars_new, JuMP.variable_by_name(pm.model, string("0_vm[",g,"]")))
            end

            load_vars = collect(variable_loads)
            for d in load_vars
                push!(vars_new, JuMP.variable_by_name(pm.model, string("0_pd[",d,"]")))
            end

            x_hat = directed_walk_ops[counter]
            N = length(vars_new)

            # additional seperating hyperplanes constraints
            @variable(pm.model, r); # r is the radius of the hypersphere from the sampled point to the relaxation
            @variable(pm.model, aux_new[1:N]) # add auxiliary variable (x_opt)
            @constraint(pm.model, aux_new .== vars_new) # constrain auxiliary variables
            @constraint(pm.model, con_sphere, sqrt(sum((aux_new[i]-x_hat[i])^2 for i in 1:N)) <= r) # add constraint on radius

            @objective(pm.model, Min, r); # the objective is to minimize r

            # check ac correction feasibility
            PF_res2 = optimize_model!(pm, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
            acpfcorrect_feasibility = nothing
            vm_vio_over = 0.0
            vm_vio_under = 0.0
            sm_vio = 0.0
            pg_vio = 0.0
            qg_vio = 0.0

            # op placeholder
            op_flag_correct = Dict(
                "N0" => 1, # flag for base-case feasibility
                "N0P" => 0.0, # active power violation
                "N0Q" => 0.0, # active power violation
                "N0OV" => 0.0, # amount over voltage violation
                "N0UV" => 0.0, # amount under voltage violation
                "N0L" => 0.0, # amount line flow violation
                "N1" => 1, # flag for N1 feasibility
                "N1OV" => 0.0, # amount over voltage violation
                "N1UV" => 0.0, # amount under voltage violation
                "N1L" => 0.0 # amount line flow violation
            )


            for i in 0:(length(multinetwork["nw"])-1)
                acpfcorrect_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res2["solution"]["nw"]["$i"], tollerance)

                if i == 0
                    if acpfcorrect_feasibility != true # || PF_res2["termination_status"] == LOCALLY_SOLVED == true  
                        op_flag_correct["N0"] = 0
                        op_flag_correct["N0P"] += pg_vio
                        op_flag_correct["N0Q"] += qg_vio
                        op_flag_correct["N0OV"] += vm_vio_over
                        op_flag_correct["N0UV"] += vm_vio_under
                        op_flag_correct["N0L"] += sm_vio
                    end
                end

                if i != 0
                    if acpfcorrect_feasibility != true # || PF_res2["termination_status"] == LOCALLY_SOLVED != true  
                        op_flag_correct["N1"] = 0
                        op_flag_correct["N1OV"] += vm_vio_over
                        op_flag_correct["N1UV"] += vm_vio_under
                        op_flag_correct["N1L"] += sm_vio
                    end
                end
            end

            # check small signal stability of perturbered generator and obtain stability indices
            sys_studied = create_system(data_opf_verif)
            construct_dynamic_model(sys_studied, dir_dynamics, case_name)
            stability = small_signal_module(sys_studied)

            if op_flag_correct["N0"] == 1 && op_flag_correct["N1"] == 1 #&& 0.0275 <= stability["damping"] <= 0.0325
                # add feasible HIC sample
                nb_feasible_dws += 1 
                push!(pf_results_feas, PF_res2["solution"]["nw"]["0"])
                push!(load_results_feas, data_opf_verif)
                push!(op_info_feas, op_flag_correct)
                push!(op_info_feas_stability, (stability["damping"], stability["distance"]))
                print("acpfcorrect feasibility: ", acpfcorrect_feasibility, "\n")
                println("ACPF correct status:", PF_res2["termination_status"] , "\n")
                print("corrected HIC point stability: ", stability["damping"], "\n")

                # add infeasible HIC sample
                push!(pf_results_infeas, PF_res0["solution"])
                push!(load_results_infeas, data_opf_verif)
                push!(op_info_infeas, op_flag)
                push!(op_info_infeas_stability, directed_walk_stability[counter])
                nb_infeasible_dws += 1
                println("initial status:", PF_res0["termination_status"] , "\n")
                print("initial feasibility: ", initial_feasibility, "\n")
                print("original HIC point stability: ", directed_walk_stability[counter][1], "\n")
            end
        end
    end

    print("number of feasible DW samples: ", nb_feasible_dws, "\n")

    # get list of feasible operating points x
    feasible_ops_dws = []
    damp_dws_feas = []
    dist_dws_feas = []
    for i in 1:length(pf_results_feas[:])
        op = get_Full_OP(data_tight_tmp, pf_results_feas[i], load_results_feas[i])
        push!(feasible_ops_dws, op)
        push!(damp_dws_feas, op_info_feas_stability[(i)][1])
        push!(dist_dws_feas, op_info_feas_stability[(i)][2])
    end

    # get list of infeasible operating points
    infeasible_ops_dws = []
    damp_dws_infeas = []
    dist_dws_infeas = []
    for i in 1:length(pf_results_infeas[:])
        op = get_Full_OP(data_tight_tmp, pf_results_infeas[i], load_results_infeas[i])
        push!(infeasible_ops_dws, op)
        push!(damp_dws_infeas, op_info_infeas_stability[(i)][1])
        push!(dist_dws_infeas, op_info_infeas_stability[(i)][2])
    end


    return feasible_ops_dws, pf_results_feas, op_info_feas, infeasible_ops_dws, pf_results_infeas, op_info_infeas, nb_feasible_dws, nb_infeasible_dws, initial_feasible_dws, damp_dws_feas, damp_dws_infeas, dist_dws_feas, dist_dws_infeas


end












# function dw_ops_feasibility(network_basic, data_tight_tmp, variable_loads, directed_walk_ops, directed_walk_stability)

#     pm, N, vars, header = instantiate_system_QCRM(data_tight_tmp, variable_loads)
#     pg_numbers, vm_numbers, pd_numbers = extract_number_and_type(vcat(header[1]))

#     pf_results_feas = []
#     load_results_feas = []
#     op_info_feas = []
#     op_info_feas_stability = []

#     pf_results_infeas = []
#     load_results_infeas = []
#     op_info_infeas = []
#     op_info_infeas_stability = []
    
#     nb_feasible_dws = 0
#     nb_infeasible_dws = 0
#     initial_feasible_dws = 0
#     tollerance = 1e-4

#     # avoid modifying original data
#     data_opf_verif = deepcopy(data_tight_tmp)

#     for counter in 1:(length(directed_walk_ops[:])) 
        
#         for g in eachindex(pg_numbers)
#             data_opf_verif["gen"]["$(pg_numbers[g])"]["pg"] = directed_walk_ops[counter][g]
#         end
#         for v in eachindex(vm_numbers)
#             data_opf_verif["bus"]["$(vm_numbers[v])"]["vm"] = directed_walk_ops[counter][length(pg_numbers)+v] 
#         end
#         for d in eachindex(pd_numbers)
#             var_load_index = pd_numbers[d]
#             data_opf_verif["load"]["$(pd_numbers[d])"]["pd"] = directed_walk_ops[counter][length(pg_numbers)+length(vm_numbers)+var_load_index] 
#             pf = data_opf_verif["load"]["$(pd_numbers[d])"]["pf"]
#             pd = data_opf_verif["load"]["$(pd_numbers[d])"]["pd"]
#             sd = pd / pf
#             qd = sqrt(sd^2 - pd^2)
#             data_opf_verif["load"]["$(pd_numbers[d])"]["qd"] = qd
#         end

#         print("Current DW sample number: ", Int(counter), "\n")

#         # dictionary placeholder with OP flags. 1 is feasible, 0 is infeasible
#         op_flag = Dict(
#             "N0" => 1, # flag for base-case feasibility
#             "N0P" => 0.0, # active power violation
#             "N0Q" => 0.0, # active power violation
#             "N0OV" => 0.0, # amount over voltage violation
#             "N0UV" => 0.0, # amount under voltage violation
#             "N0L" => 0.0, # amount line flow violation
#             "N1" => 1, # flag for N1 feasibility
#             "N1OV" => 0.0, # amount over voltage violation
#             "N1UV" => 0.0, # amount under voltage violation
#             "N1L" => 0.0 # amount line flow violation
#         )

#         # construct SCOPF in form of multi network formulation
#         multinetwork = build_c1_scopf_multinetwork_modif(data_opf_verif)
#         # if length(multinetwork["nw"]) == 1
#         #     op_flag["N1"] = 0
#         # end

#         if multinetwork["per_unit"] == true
#             for (n, network) in multinetwork["nw"]
#                 multinetwork["nw"]["$n"]["per_unit"] = true
#             end
#         end

#         # check initial feasibility of base case and contingency cases
#         PF_res0 = nothing
#         initial_feasibility = nothing
#         vm_vio_over = 0.0
#         vm_vio_under = 0.0
#         sm_vio = 0.0

#         for i in 0:(length(multinetwork["nw"])-1)
#             PF_res0 = solve_ac_pf(multinetwork["nw"]["$i"], optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
#             update_data!(multinetwork["nw"]["$i"], PF_res0["solution"]) # update data with PF results
#             flows0 = calc_branch_flow_ac(multinetwork["nw"]["$i"]) # compute branch flows
#             update_data!(multinetwork["nw"]["$i"], flows0) # add branch flows
#             update_data!(PF_res0["solution"], flows0) # add branch flows to solution
#             initial_feasibility, pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio = check_ac_feasibility(data_tight_tmp, PF_res0["solution"], tollerance)

#             print(pg_vio, qg_vio, vm_vio_over, vm_vio_under, sm_vio)
#             if i == 0
#                 # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # initial_feasibility != true
#                 if initial_feasibility != true 
#                     op_flag["N0"] = 0
#                     op_flag["N0P"] += pg_vio
#                     op_flag["N0Q"] += qg_vio
#                     op_flag["N0OV"] += vm_vio_over
#                     op_flag["N0UV"] += vm_vio_under
#                     op_flag["N0L"] += sm_vio
#                 end
#             end

#             if i != 0
#                 # if PF_res0["termination_status"] != LOCALLY_SOLVED  || PF_res0["primal_status"] != FEASIBLE_POINT || PF_res0["dual_status"] != FEASIBLE_POINT # initial_feasibility != true
#                 if initial_feasibility != true 
#                     op_flag["N1"] = 0
#                     op_flag["N1OV"] += vm_vio_over
#                     op_flag["N1UV"] += vm_vio_under
#                     op_flag["N1L"] += sm_vio
#                 end
#             end

#         end                

#         # check feasibility
#         if op_flag["N0"] == 1 && op_flag["N1"] == 1
#             nb_feasible_dws += 1
#             initial_feasible_dws += 1
#             push!(pf_results_feas, PF_res0["solution"])
#             push!(load_results_feas, data_opf_verif)
#             push!(op_info_feas, op_flag) 
#             push!(op_info_feas_stability, directed_walk_stability[counter])
#             println("initial status:", PF_res0["termination_status"] , "\n")
#             print("initial feasibility: ", initial_feasibility, "\n")
#         else 
#             # add infeasible initial sample
#             push!(pf_results_infeas, PF_res0["solution"])
#             push!(load_results_infeas, data_opf_verif)
#             push!(op_info_infeas, op_flag)
#             push!(op_info_infeas_stability, directed_walk_stability[counter])
#             nb_infeasible_dws += 1
#             println("initial status:", PF_res0["termination_status"] , "\n")
#             print("initial feasibility: ", initial_feasibility, "\n")
#         end
#     end

#     print("number of feasible DW samples: ", nb_feasible_dws, "\n")

#     # get list of feasible operating points x
#     feasible_ops_dws = []
#     damp_dws_feas = []
#     dist_dws_feas = []
#     for i in 1:length(pf_results_feas[:])
#         op = get_Full_OP(data_tight_tmp, pf_results_feas[i], load_results_feas[i])
#         push!(feasible_ops_dws, op)
#         push!(damp_dws_feas, op_info_feas_stability[(i)][1])
#         push!(dist_dws_feas, op_info_feas_stability[(i)][2])
#     end

#     # get list of infeasible operating points
#     infeasible_ops_dws = []
#     damp_dws_infeas = []
#     dist_dws_infeas = []
#     for i in 1:length(pf_results_infeas[:])
#         op = get_Full_OP(data_tight_tmp, pf_results_infeas[i], load_results_infeas[i])
#         push!(infeasible_ops_dws, op)
#         push!(damp_dws_infeas, op_info_infeas_stability[(i)][1])
#         push!(dist_dws_infeas, op_info_infeas_stability[(i)][2])
#     end


#     return feasible_ops_dws, pf_results_feas, op_info_feas, infeasible_ops_dws, pf_results_infeas, op_info_infeas, nb_feasible_dws, nb_infeasible_dws, initial_feasible_dws, damp_dws_feas, damp_dws_infeas, dist_dws_feas, dist_dws_infeas


# end






















