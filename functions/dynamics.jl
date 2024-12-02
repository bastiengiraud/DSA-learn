using CSV
using DataFrames
using PowerModels
using PowerSystems
using UUIDs

# Function to convert string to number or keep it as is
function convert_value(value)
    if typeof(value) == String7
        value = string(value)  # Convert String7 to String
    end

    # Try to convert to Float64
    try
        return parse(Float64, value)
    catch e
        return value  # If conversion fails, return the original value
    end
end

function df_to_dict(df)
    # Extract the column names for the initial keys
    initial_keys = df[1, 2:end]

    # Create the nested dictionary
    nested_dict = Dict{String, Dict{String, Any}}()

    # Populate the nested dictionary
    for col in 2:ncol(df)
        unit_number = df[1, col]
        nested_dict[string(unit_number)] = Dict{String, Any}()
        for row in 2:nrow(df)
            key = df[row, 1]
            value = df[row, col]
            nested_dict[string(unit_number)][key] = convert_value(value)
        end
    end
    
    return nested_dict
end

function load_machine_data(dir_dynamics, case_name)
    # Define file paths
    dir_machine = joinpath(dir_dynamics, case_name, "machine.csv")
    dir_avr = joinpath(dir_dynamics, case_name, "avr.csv")
    dir_gov = joinpath(dir_dynamics, case_name, "gov.csv")
    dir_pss = joinpath(dir_dynamics, case_name, "pss.csv")

    # Read data and convert to dictionaries
    machine_data = df_to_dict(CSV.read(dir_machine, DataFrame; delim=';', header=false))
    avr_data = df_to_dict(CSV.read(dir_avr, DataFrame; delim=';', header=false))
    gov_data = df_to_dict(CSV.read(dir_gov, DataFrame; delim=';', header=false))
    pss_data = df_to_dict(CSV.read(dir_pss, DataFrame; delim=';', header=false))

    # Combine all into a single dictionary
    machine_data_dict = Dict(
        :machine => machine_data,
        :avr => avr_data,
        :gov => gov_data,
        :pss => pss_data
    )

    return machine_data_dict
end






function create_machine(machine_data, gen_number)
    # https://nrel-sienna.github.io/PowerSystems.jl/v1.8/model_library/generated_Machine/#RoundRotorQuadratic
    
    operation = machine_data["$gen_number"]["Operation"]

    if machine_data["$gen_number"]["Type"] == "GENROU"

        return machine_genrou = RoundRotorQuadratic(
        R = machine_data["$gen_number"]["ra"],
        Td0_p = machine_data["$gen_number"]["Td0_p"],
        Td0_pp = machine_data["$gen_number"]["Td0_pp"],
        Tq0_p = machine_data["$gen_number"]["Tq0_p"],
        Tq0_pp = machine_data["$gen_number"]["Tq0_pp"],
        Xd = machine_data["$gen_number"]["Xd"],
        Xq = machine_data["$gen_number"]["Xq"],
        Xd_p = machine_data["$gen_number"]["Xd_p"],
        Xq_p = machine_data["$gen_number"]["Xq_p"],
        Xd_pp = machine_data["$gen_number"]["Xd_pp"], # Xd_pp = Xq_pp
        Xl = machine_data["$gen_number"]["Xl"],
        Se = (machine_data["$gen_number"]["S1"], machine_data["$gen_number"]["S1.2"]) 
        ), operation

    elseif machine_data["$gen_number"]["Type"] == "ODOQ"

            return machine_odoq = OneDOneQMachine(
            R = machine_data["$gen_number"]["ra"],
            Td0_p = machine_data["$gen_number"]["Td0_p"],
            Tq0_p = machine_data["$gen_number"]["Tq0_p"],
            Xd = machine_data["$gen_number"]["Xd"],
            Xq = machine_data["$gen_number"]["Xq"],
            Xd_p = machine_data["$gen_number"]["Xd_p"],
            Xq_p = machine_data["$gen_number"]["Xq_p"],
            ), operation

    end


end

function create_shaft(machine_data, gen_number) # 39-bus machine old


    if machine_data["$gen_number"]["Type"] == "GENROU"

        return shaft_no_damping_genrou = SingleMass(
            H = machine_data["$gen_number"]["H"],
            D = machine_data["$gen_number"]["D"],
        )

    elseif machine_data["$gen_number"]["Type"] == "ODOQ"

        return shaft_no_damping_genrou = SingleMass(
            H = machine_data["$gen_number"]["H"],
            D = machine_data["$gen_number"]["D"],
        )

    end

end

#avrdata = create_avr(dir_dynamics, case_name, 1)

function create_avr(avr_data, gen_number) # 39-bus avr old

    if avr_data["$gen_number"]["Type"] == "AVRTypeI" # 

        #return avr_none = AVRFixed(0.0)
        return avr = PowerSystems.AVRTypeI(
            avr_data["$gen_number"]["Ka"], #Ka - Gain
            avr_data["$gen_number"]["Ke"], #Ke
            avr_data["$gen_number"]["Kf"], #Kf
            avr_data["$gen_number"]["Ta"], #Ta
            avr_data["$gen_number"]["Te"], #Te
            avr_data["$gen_number"]["Tf"], #Tf
            avr_data["$gen_number"]["Tr"], #Tr
            (min = avr_data["$gen_number"]["Vmin"], max = avr_data["$gen_number"]["Vmax"]),
            avr_data["$gen_number"]["Ae"], #Ae - 1st ceiling coefficient
            avr_data["$gen_number"]["Be"], #Be - 2nd ceiling coefficient
            )

    elseif avr_data["$gen_number"]["Type"] == "AVRTypeII" # 

            # return avr_none = AVRFixed(0.0)
            return avr = PowerSystems.AVRTypeII(
                avr_data["$gen_number"]["Ka"], #Ka - Gain
                avr_data["$gen_number"]["Ke"], #Ke
                avr_data["$gen_number"]["Kf"], #Kf
                avr_data["$gen_number"]["Ta"], #Ta
                avr_data["$gen_number"]["Te"], #Te
                avr_data["$gen_number"]["Tf"], #Tf
                avr_data["$gen_number"]["Tr"], #Tr
                (min = avr_data["$gen_number"]["Vmin"], max = avr_data["$gen_number"]["Vmax"]),
                avr_data["$gen_number"]["Ae"], #Ae - 1st ceiling coefficient
                avr_data["$gen_number"]["Be"], #Be - 2nd ceiling coefficient
                )

    elseif avr_data["$gen_number"]["Type"] == "EXST1"

        #return avr_none = AVRFixed(0.0)
        return avr_exst1 = PowerSystems.EXST1(
            Tr = avr_data["$gen_number"]["Tr"],
            Vi_lim=(min=avr_data["$gen_number"]["Vimin"], max=avr_data["$gen_number"]["Vimax"]),
            Tc = avr_data["$gen_number"]["Tc"],
            Tb = avr_data["$gen_number"]["Tb"],
            Ka = avr_data["$gen_number"]["Ka"],
            Ta = avr_data["$gen_number"]["Ta"],
            Vr_lim=(min=avr_data["$gen_number"]["Vrmin"], max=avr_data["$gen_number"]["Vrmax"]),
            Kc = avr_data["$gen_number"]["Kc"],
            Kf = avr_data["$gen_number"]["Kf"],
            Tf = avr_data["$gen_number"]["Tf"],
            #V_ref = avr_data["$gen_number"]["Vref"],
            #ext=Dict{String, Any}(),
            )

    elseif avr_data["$gen_number"]["Type"] == "NONE"

        return avr_none = AVRFixed(0.0)

    else
        return avr_none = AVRFixed(0.0)

    end

end


function create_gov(gov_data, operation, gen_number) # 39-bus: gov 2

    if operation == "NG"

        if gov_data["$gen_number"]["Type"] == "TGTypeI" 

            #return tg_none = TGFixed(efficiency = 1.0) 
            return tg_type1 = TGTypeI(
            gov_data["$gen_number"]["R"], #R, droop parameter
            gov_data["$gen_number"]["Ts"], #Ts, governor time constant
            gov_data["$gen_number"]["Tc"], #Tc, servo time constant
            gov_data["$gen_number"]["T3"], #T3, transient gain time constant
            gov_data["$gen_number"]["T4"], #T4, power fraction time constant
            gov_data["$gen_number"]["T5"], #T5, reheat time constant
            (min = gov_data["$gen_number"]["min"], max = gov_data["$gen_number"]["max"]), #P_lims
            )

        elseif gov_data["$gen_number"]["Type"] == "TGOVI" 

            #return tg_none = TGFixed(efficiency = 1.0) 
            return tg_type1 = SteamTurbineGov1(
            R = gov_data["$gen_number"]["R"], #R, droop parameter
            T1 = gov_data["$gen_number"]["T1"], #Ts, governor time constant
            T2 = gov_data["$gen_number"]["T2"], #Tc, servo time constant
            T3 = gov_data["$gen_number"]["T3"], #T3, transient gain time constant
            D_T = gov_data["$gen_number"]["DT"], #T4, power fraction time constant
            DB_h = 0,
            DB_l = 0,
            T_rate = 0,
            valve_position_limits = (min = gov_data["$gen_number"]["Vmin"], max = gov_data["$gen_number"]["Vmax"]), #P_lims
            )   

        else

            return tg_none = TGFixed(1.0) 

        end

    elseif operation == "SYNC"
        return tg_none = TGFixed(1.0) 

    elseif operation == "CON"
        return tg_none = TGFixed(1.0) 

    else
        return tg_none = TGFixed(1.0) 
    end

end


function create_pss(pss_data, operation, gen_number)
    # http://www1.sel.eesc.usp.br/ieee/TwoArea/Two%20area%20four%20generator%20model%20PSSE%20study%20report.pdf

    if pss_data["$gen_number"]["Type"] == "IEEEST" 

        #return pss_none = PSSFixed(V_pss = 0.0)
        return pss = IEEEST(
            input_code = pss_data["$gen_number"]["input_code"],
            remote_bus_control = pss_data["$gen_number"]["remote_bus_control"],
            A1 = pss_data["$gen_number"]["A1"],
            A2 = pss_data["$gen_number"]["A2"],
            A3 = pss_data["$gen_number"]["A3"],
            A4 = pss_data["$gen_number"]["A4"],
            A5 = pss_data["$gen_number"]["A5"],
            A6 = pss_data["$gen_number"]["A6"],
            T1 = pss_data["$gen_number"]["T1"],
            T2 = pss_data["$gen_number"]["T2"],
            T3 = pss_data["$gen_number"]["T3"],
            T4 = pss_data["$gen_number"]["T4"],
            T5 = pss_data["$gen_number"]["T5"],
            T6 = pss_data["$gen_number"]["T6"],
            Ks = pss_data["$gen_number"]["Ks"],
            Ls_lim = (pss_data["$gen_number"]["Lsmin"], pss_data["$gen_number"]["Lsmax"]),
            Vcu = pss_data["$gen_number"]["Vcu"],
            Vcl = pss_data["$gen_number"]["Vcl"],
            #ext=Dict{String, Any}(),
            )

    else

        return pss_none = PSSFixed(0.0)
    end
end


function initialize_steady_state_model(file_path)

    system_model = System(file_path)

    return system_model
end


# function create_system(network)
#     export_matpower("file.m", network)
#     file_path = "file.m"
#     # The target line to detect and comment out
#     target_line = "mpc.multinetwork"
#     # Call the function to process the file
#     comment_line_in_file(file_path, target_line)
#     target_line = "mpc.multiinfrastructure"
#     # Call the function to process the file
#     comment_line_in_file(file_path, target_line)

#     sys = System("file.m") # function from PowerSystems

#     return sys
# end

function create_system(network)
    file_id = UUIDs.uuid4()  # Generate a random UUID
    file_path = "file_$file_id.m"  # Use the UUID in the filename

    export_matpower(file_path, network)  # Export to the unique filename

    # Comment out target lines
    target_lines = ["mpc.multinetwork", "mpc.multiinfrastructure"]
    for target_line in target_lines
        comment_line_in_file(file_path, target_line)
    end

    sys = System(file_path)  # Use the unique file path

    rm(file_path) # remove the file after usage

    return sys
end



function construct_dynamic_model(system, machine_data_dict)

    # Access dictionary elements using keys
    machine_data = machine_data_dict[:machine]
    avr_data = machine_data_dict[:avr]
    gov_data = machine_data_dict[:gov]
    pss_data = machine_data_dict[:pss]

    # add NGs and SYNCs 
    for g in get_components(Generator, system)
        # get generator bus
        generator_bus = get_number(get_bus(g))

        # make dynamic data for the generator at that bus
        machine, operation = create_machine(machine_data, generator_bus)
        shaft = create_shaft(machine_data, generator_bus)
        avr = create_avr(avr_data, generator_bus)
        gov = create_gov(gov_data, operation, generator_bus)
        pss = create_pss(pss_data, operation, generator_bus)

        # construct dynamic generator
        case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine,
            shaft = shaft,
            avr = avr,
            prime_mover = gov,
            pss = pss,
        )

        #Attach the dynamic generator to the system
        add_component!(system, case_gen, g)

    end
end


#system_created = create_system(Initialize.network_data)
#dynsys = construct_dynamic_model(system_created, Initialize.dir_dynamics, Initialize.case_name)
#machine, operation, machine_data = create_machine(dir_dynamics, case_name, 16)

# to_json(threebus_sys, "YOUR_DIR/threebus_sys.json")




















