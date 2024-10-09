
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
        if isempty(zero_real_and_imag_indices)
            println("Warning : 0 but stable")
        else
            println("Warning : there is an eigenvalue at the origin")
        end
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
    ss_stability = false
    sys = full_model
    time_span = (0.0, 30.0)
    sim_studied = Simulation(ResidualModel, sys , pwd(), time_span) # , initialize_simulation = true)#, initial_conditions = initial_conditions) # origineel
    # incon = read_initial_conditions(sim_studied)
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

    return stability
end


function dynamic_components()

    #Dynamic components creation
    machine_oneDoneQ() = OneDOneQMachine(
        0.0, #R
        0.146, #Xd
        0.0969, #Xq
        0.0608, #Xd_p
        0.0969, #Xq_p
        8.96, #Td0_p
        0.31, #Tq0_p
    )

    machine_twoDoneQ() = OneDOneQMachine(
        0.0, #R
        0.8958, #Xd
        0.8645, #Xq
        0.1198, #Xd_p
        0.1969, #Xq_p
        6.0, #Td0_p
        0.535, #Tq0_p
    )

    machine_threeDoneQ() = OneDOneQMachine(
        0.0, #R
        1.3125, #Xd
        1.2578, #Xq
        0.1813, #Xd_p
        0.25, #Xq_p
        5.89, #Td0_p
        0.6, #Tq0_p
    )

    machine_keys = [1, 2, 3]
    machine_values = [machine_oneDoneQ(), machine_twoDoneQ(), machine_threeDoneQ()]
    machine = Dict(zip(machine_keys, machine_values))

    #shaft_special() = FiveMassShaft(23.64,0,0,0,0,0.125,0,0,0,0,0,0,0,0,0,0,0,0)
    shaft_one() = SingleMass(
        23.64, #H 
        0.125, #D
    )

    shaft_two() = SingleMass(
        6.4, #H 
        0.033, #D
    )

    shaft_three() = SingleMass(
        3.01, #H 
        0.016, #D
    )

    shaft_keys = [1, 2, 3]
    shaft_values = [shaft_one(), shaft_two(), shaft_three()]
    shaft = Dict(zip(shaft_keys, shaft_values))

    avr_type1() = AVRTypeI(
        20.0, #Ka - Gain
        1.0, #Ke
        0.063, #Kf
        0.2, #Ta
        0.314, #Te
        0.35, #Tf
        0.001, #Tr
        (min = -5.0, max = 5.0),
        0.0039, #Ae - 1st ceiling coefficient
        1.555, #Be - 2nd ceiling coefficient
    )

    #No TG
    tg_none() = TGFixed(1.0) #efficiency

    #No PSS
    pss_none() = PSSFixed(0.0) #Vs

    generator_keys = ["machine", "shaft", "AVR", "TG", "PSS"]
    generator_values = [machine, shaft, avr_type1(), tg_none(), pss_none()]
    generator = Dict(zip(generator_keys, generator_values))

    # inverter
    converter_high_power = AverageConverter(rated_voltage = 13.80, rated_current = 100.0)
    outer_cont = OuterControl(
            VirtualInertia(Ta = 2.0, kd = 400.0, kω = 20.0),
            ReactivePowerDroop(kq = 0.2, ωf = 1000.0),
        )
    inner_cont = VoltageModeControl(
            kpv = 0.59,     #Voltage controller proportional gain
            kiv = 736.0,    #Voltage controller integral gain
            kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
            rv = 0.0,       #Virtual resistance in pu
            lv = 0.2,       #Virtual inductance in pu
            kpc = 1.27,     #Current controller proportional gain
            kic = 14.3,     #Current controller integral gain
            kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
            ωad = 50.0,     #Active damping low pass filter cut-off frequency
            kad = 0.2,      #Active damping gain
        )   
    dc_source_lv = FixedDCSource(voltage = 600.0)   

    pll = KauraPLL(
            ω_lp = 500.0, #Cut-off frequency for LowPass filter of PLL filter.
            kp_pll = 0.084,  #PLL proportional gain
            ki_pll = 4.69,   #PLL integral gain
        )
    filt = LCLFilter(lf = 0.08, rf = 0.003, cf = 0.074, lg = 0.2, rg = 0.01)

    converter_keys = ["converter", "outer_control", "inner_control", "dc_source", "PLL", "filt"]
    converter_values = [converter_high_power, outer_cont, inner_cont, dc_source_lv, pll, filt]
    converter = Dict(zip(converter_keys, converter_values))

    # create components dict
    dynamic_components_dict_keys = ["generator", "converter"]
    dynamic_components_dict_values = [generator, converter]
    dynamic_components_dict = Dict(zip(dynamic_components_dict_keys, dynamic_components_dict_values))

    return dynamic_components_dict 
end


function initialize_steady_state_model(file_path)

    system_model = System(file_path)

    return system_model
end



function construct_model(system_model, gen_bus)

    slack_bus = [b for b in get_components(Bus, system_model) if b.bustype == BusTypes.REF][1]

    dynamic_components_dict = dynamic_components()

    counter = 0
    previous_gen = 0
    for i in gen_bus
        counter += 1
        for g in get_components(Generator, system_model) # get number of generators in system
            # print(g) # g is all the info of all generators; ThermalStandard, LoadZone, Bus, ThreePartCost, PrimeMovers, ThermalFuels
            # println(get_bus(g)) # get_bus gets all the bus info; bus number, REF, voltage, limits area etc etc
            # println(get_number(get_bus(g))) # get_number gets the bus number

            bus_number = get_number(get_bus(g))
            if bus_number == gen_bus[counter] && get_name(g) != previous_gen # the bus number should match the bus index
                # also make sure the for loop skips any generators which are already connected to the bus
                # gen_bus should be in ascending order, so buses with multiple generators should skip the previous
                #Create the dynamic generator
                gen_type = rand(1:3) # for now, assign a random gen type to each generator just for demo purposes
                case_gen = DynamicGenerator(
                    get_name(g),
                    1.0, # ω_ref,
                    dynamic_components_dict["generator"]["machine"][gen_type], # machine_oneDoneQ(), #machine
                    dynamic_components_dict["generator"]["shaft"][gen_type], # shaft_one(), #shaft
                    dynamic_components_dict["generator"]["AVR"], # avr_type1(), #avr
                    dynamic_components_dict["generator"]["TG"], # tg_none(), #tg
                    dynamic_components_dict["generator"]["PSS"], # , #pss
                )
                #Attach the dynamic generator to the system by specifying the dynamic and static part
                add_component!(system_model, case_gen, g)

                # ensure you skip generator if multiple are connected to the same bus
                gen_bus[counter] = 0
                previous_gen = get_name(g)
                
            end

        end
    end

    return system_model
end









##---------- below here moved from module_methods


machine_oneDoneQ() = OneDOneQMachine(
           R = 0.0,
           Xd = 1.3125,
           Xq = 1.2578,
           Xd_p = 0.1813,
           Xq_p = 0.25,
           Td0_p = 5.89,
           Tq0_p = 0.6,
       )

machine1() = RoundRotorQuadratic(
R = 0.0043 ,
Td0_p = 0.5871,
Td0_pp = 0.0248,
Tq0_p = 0.1351,
Tq0_pp = 0.0267,
Xd = 1.670,
Xq = 1.600,
Xd_p = 0.265,
Xq_p = 0.460,
Xd_pp = 0.205,
Xl = 0.150,
Se = (0.091, 0.400 )
)

avr_type1() = AVRTypeI(
           Ka = 20.0,
           Ke = 0.01,
           Kf = 0.063,
           Ta = 0.2,
           Te = 0.314,
           Tf = 0.35,
           Tr = 0.001,
           Va_lim = (min = -5.0, max = 5.0),
           Ae = 0.0039, #1st ceiling coefficient
           Be = 1.555, #2nd ceiling coefficient
       )

tg_none() = TGFixed(efficiency = 1.0)       
pss_none() = PSSFixed(V_pss = 0.0)
shaft_no_damping() = SingleMass(
           H = 3.01,
           D = 0.0,
       )

TGOV1() = SteamTurbineGov1(
    R = 0.05,
    T1 = 0.49,
    valve_position_limits = (min = 33 , max = 0.4),
    T2 = 2.1,
    T3 = 7,
    D_T = 0,
    DB_h = 0,
    DB_l = 0,
    T_rate = 0,
)       


tg_type1() = TGTypeI(
0.02, #R
0.1, #Ts
0.45, #Tc
0.0, #T3
12.0, #T4
50.0, #T5
(min = 0.0, max = 1.2), #P_lims
)

machine1() = RoundRotorQuadratic(
    R = 0.0043 ,
    Td0_p = 0.5871,
    Td0_pp = 0.0248,
    Tq0_p = 0.1351,
    Tq0_pp = 0.0267,
    Xd = 1.670,
    Xq = 1.600,
    Xd_p = 0.265,
    Xq_p = 0.460,
    Xd_pp = 0.205,
    Xl = 0.150,
    Se = (0.091, 0.400 )
)

machine2() = RoundRotorQuadratic(
    R = 0.0035 ,
    Td0_p = 1.100,
    Td0_pp = 0.0277,
    Tq0_p = 0.1086 ,
    Tq0_pp = 0.0351,
    Xd = 1.180 ,
    Xq = 1.050 ,
    Xd_p = 0.220,
    Xq_p = 0.380 ,
    Xd_pp = 0.145 ,
    Xl = 0.075 ,
    Se = (0.0933, 0.4044  )
)

machine3() = RoundRotorQuadratic(
    R = 0.0 ,
    Td0_p = 11.600 ,
    Td0_pp = 0.058,
    Tq0_p = 0.159 ,
    Tq0_pp = 0.201,
    Xd = 2.373,
    Xq = 1.172 ,
    Xd_p = 0.343,
    Xq_p = 1.172,
    Xd_pp = 0.231 ,
    Xl = 0.150,
    Se = (0.091, 0.400 )
)

machine4() = RoundRotorQuadratic(
    R = 0.0025 ,
    Td0_p = 8.000,
    Td0_pp = 0.0525,
    Tq0_p = 0.008 ,
    Tq0_pp = 0.0151,
    Xd = 1.769 ,
    Xq = 0.855 ,
    Xd_p = 0.304 ,
    Xq_p = 0.5795,
    Xd_pp = 0.2035 ,
    Xl = 0.1045,
    Se = (0.304, 0.666)
)

avr_none() = AVRFixed(0.0)

machine_genrou() = RoundRotorExponential(;
    R = 0.0,
    Td0_p = 8.0,
    Td0_pp = 0.03,
    Tq0_p = 0.4,
    Tq0_pp = 0.05,
    Xd = 1.8,
    Xq = 1.7,
    Xd_p = 0.3,
    Xq_p = 0.55,
    Xd_pp = 0.25,
    Xl = 0.2,
    Se = (0.0, 0.0),
)

shaft_no_damping_genrou() = SingleMass(
           H = 6.5,
           D = 0.0,
       )

avr_type1() = AVRTypeI(
           Ka = 20.0,
           Ke = 0.01,
           Kf = 0.063,
           Ta = 0.2,
           Te = 0.314,
           Tf = 0.35,
           Tr = 0.001,
           Va_lim = (min = -5.0, max = 5.0),
           Ae = 0.0039, #1st ceiling coefficient
           Be = 1.555, #2nd ceiling coefficient
       )

sexs() = SEXS(
    Ta_Tb = 0.1,
    Tb = 10,
    K = 100,
    Te = 0.1,
    V_lim = (min = 0.0, max = 5.0),
)

tg_none() = TGFixed(efficiency = 1.0)       
pss_none() = PSSFixed(V_pss = 0.0)

shaft1() = SingleMass(H = 2.656 , D = 2.0)
shaft2() = SingleMass(H = 4.985 , D = 2.0)
shaft3() = SingleMass(H = 1.520 , D = 0.0)
shaft4() = SingleMass(H = 1.20 , D = 0.0)

tg_type1() = TGTypeI(
    0.02, #R
    0.1, #Ts
    0.45, #Tc
    0.0, #T3
    12.0, #T4
    50.0, #T5
    (min = 0.0, max = 1.2), #P_lims
)

function construct_dyn_14bus(sys_14)
    for g in get_components(Generator, sys_14)
        if get_number(get_bus(g)) == 1
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine1(),
            shaft = shaft1(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)

        end 

        if get_number(get_bus(g)) == 2
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine2(),
            shaft = shaft2(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)

        end 

        if get_number(get_bus(g)) == 3
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine3(),
            shaft = shaft3(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)

        end 

        if get_number(get_bus(g)) == 6  || get_number(get_bus(g)) == 8
            case_gen = DynamicGenerator(
            name = get_name(g),
            ω_ref = 1.0,
            machine = machine4(),
            shaft = shaft4(),
            avr = avr_type1(),
            prime_mover = tg_type1(),
            pss = pss_none(),
        )
            #Attach the dynamic generator to the system
            add_component!(sys_14, case_gen, g)
        end 
    end
end




