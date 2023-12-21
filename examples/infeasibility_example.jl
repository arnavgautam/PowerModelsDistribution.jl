using Pkg
Pkg.activate(pwd())
using Revise
using PowerModelsDistribution
using Ipopt
using BARON
using Gurobi
using JuMP
using BilevelJuMP
using Infiltrator
import MultiObjectiveAlgorithms as MOA
import MathOptInterface as MOI
using DelimitedFiles
using DataFrames
using CSV

#################################################################################################
# Setting up information on the case study                                                      #
#################################################################################################

const SFO_SMART_DS_REGION = "P4U"
const FILENAME_BASE = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_$(SFO_SMART_DS_REGION)_Timeseries/scenarios/base_timeseries/opendss/p4uhs0_4/p4uhs0_4--p4udt0/"

# sf_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_$(SFO_SMART_DS_REGION)_Timeseries/scenarios/base_timeseries/opendss/Master.dss"
sf_casefile = joinpath(FILENAME_BASE, "Master.dss")

# Setting up custom optimizers
# my_ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>200, "tol"=>100.0, "dual_inf_tol"=>100.0, "constr_viol_tol"=>1.0, "compl_inf_tol"=>0.01)
my_ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>5000, "acceptable_tol"=>100000.0, "acceptable_dual_inf_tol"=>10000.0, "acceptable_constr_viol_tol"=>1.0, "acceptable_compl_inf_tol"=>1.0, "acceptable_obj_change_tol"=>100.0)
tighter_ipopt_optimizer = optimizer_with_attributes(
    Ipopt.Optimizer,
    "max_iter"=>5000,
    "acceptable_tol"=>100000.0,
    "acceptable_dual_inf_tol"=>10000.0,
    "acceptable_constr_viol_tol"=>1.0,
    "acceptable_compl_inf_tol"=>1.0,
    "acceptable_obj_change_tol"=>1e-6
    )
my_gurobi_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "TimeLimit"=>60,#36000,
    "MipGap"=>2e-5,#0.01,
    "nonconvex"=>2) # "gurobi_options"=>"timelim 300 mipgap 0.05 nonconvex 2"
my_baron_optimizer = optimizer_with_attributes(
    BARON.Optimizer,
    "threads"=>1,#6,
    "EpsA"=>1e-12,#e6,
    "EpsR"=>1e-12,#e6,
    "DeltaTerm"=>1,
    "DeltaT"=>600,
    "AbsConFeasTol"=>1e6,
    "RelConFeasTol"=>0.1,
    "AbsIntFeasTol"=>1e6,
    "RelIntFeasTol"=>0.1,
    "FirstFeas"=>1,
    "FirstLoc"=>1,
    "MaxTime"=>36000
    )

#################################################################################################
# CORE FUNCTIONS                                                                                #
#################################################################################################

# Make a mathematical model to represent the electric grid in the case file
function create_PMD_mathematical_model(case_file; is_multinetwork=true, time_series_to_run="yearly")
    eng_model = parse_file(case_file; multinetwork=is_multinetwork, time_series=time_series_to_run)
    set_time_elapsed!(eng_model, 0.25)

    # Remove the infinite power injection of the slack bus
    # eng_model["voltage_source"]["source"]["status"]=DISABLED
    # eng_model["voltage_source"]["source"]["rs"] .= 1e100#Inf
    # eng_model["voltage_source"]["source"]["xs"] .= 1e100#Inf
    # eng_model["voltage_source"]["source"]["pg_ub"] = [0, 0, 0]
    # eng_model["voltage_source"]["source"]["qg_ub"] = [0, 0, 0]

    # for nw in 1:length(eng_model["nw"])
    #     eng_model["nw"]["$(nw-1)"]["voltage_source"]["source"]["pg_ub"] = [0, 0, 0]
    #     eng_model["nw"]["$(nw-1)"]["voltage_source"]["source"]["qg_ub"] = [0, 0, 0]
    #     eng_model["nw"]["$(nw-1)"]["voltage_source"]["source"]["pg_lb"] = [0, 0, 0]
    #     eng_model["nw"]["$(nw-1)"]["voltage_source"]["source"]["qg_lb"] = [0, 0, 0]
    #     eng_model["nw"]["$(nw-1)"]["voltage_source"]["source"]["pg"] = [0, 0, 0]
    #     eng_model["nw"]["$(nw-1)"]["voltage_source"]["source"]["qg"] = [0, 0, 0]
    #     eng_model["nw"]["$(nw-1)"]["voltage_source"]["source"]["vm"] = [Inf, Inf, Inf, Inf]
    # end

    # Convert the ENGINEERING model into a MATHEMATICAL model
    passthrough_dict =  Dict{String,Vector{String}}(
        "bus" => String["lat", "lon"])
    math_model = transform_data_model(eng_model, eng2math_passthrough=passthrough_dict)

    # for nw in 1:length(eng_model["nw"])
    #     @assert(math_model["nw"]["$(nw-1)"]["gen"]["1"]["pmax"] == [0.0,0.0,0.0])
    #     @assert(math_model["nw"]["$(nw-1)"]["gen"]["1"]["qmax"] == [0.0,0.0,0.0])
    #     @assert(math_model["nw"]["$(nw-1)"]["gen"]["1"]["vg"] == [Inf, Inf, Inf])
    # end

    return math_model
end

# Make a JuMP model of this math model using this formulation, for this power flow problem, with this optimizer
function create_PMD_JuMP_model(math_model, formulation, problem, my_optimizer; PMD_jump_model = nothing)
    if isnothing(PMD_jump_model)
        PMD_jump_model = JuMP.Model() 
    end
    instantiate_mc_model(math_model, formulation, problem; jump_model = PMD_jump_model);
    set_optimizer(PMD_jump_model, my_optimizer)
    return PMD_jump_model
end

function write_DAT_file(timeseries_data; DAT_file_name_suffix="")
    slack_power_filename = "/home/arnav.gautam/model/data/$(SFO_SMART_DS_REGION)$(DAT_file_name_suffix)_load.dat"
    slack_power_file = open(slack_power_filename, "w")
    write(slack_power_file, "param: d_p1 := \n")
    for nw_num in 1:length(timeseries_data)
        write(slack_power_file, "$nw_num  total $(timeseries_data[nw_num])\r\n")
        write(slack_power_file, "$nw_num  cl $(timeseries_data[nw_num])\r\n")
        write(slack_power_file, "$nw_num  rl 0.0\r\n")
    end
    write(slack_power_file, ";\n")
    close(slack_power_file)
    return slack_power_filename
end

# Run PMD.jl on this math model using this formulation, for this power flow problem, with this optimizer
function run_PMD_directly_and_write_result_to_file(PMD_math_model, formulation, problem, my_optimizer; is_multinetwork=true, write_to_DAT_file=true)
    sf_tpia_solution = problem(PMD_math_model, formulation, my_optimizer; multinetwork=is_multinetwork)
    num_timestamps = length(keys(sf_tpia_solution["solution"]["nw"]))
    slack_power_timeseries_data = Float64[num_timestamps]
    for nw_num in 1:num_timestamps
        solution_p_slack = sum(sum(val["p_slack_out"]) - sum(val["p_slack_in"]) for val in values(sf_tpia_solution["solution"]["nw"]["$(nw_num-1)"]["bus"]))
        solution_q_slack = sum(sum(val["q_slack_out"]) - sum(val["q_slack_in"]) for val in values(sf_tpia_solution["solution"]["nw"]["$(nw_num-1)"]["bus"]))
        solution_slack_power = (solution_p_slack^2 + solution_q_slack^2).^0.5
        slack_power_timeseries_data[nw_num] = solution_slack_power
    end
    if write_to_DAT_file
        return write_DAT_file(slack_power_timeseries_data)
    else
        return slack_power_timeseries_data
    end
end

# Run a Python script to customize an AMPL economic optimization model with load data from this region, and output it to text
function refresh_economic_model_in_AMPL(region)
    run(`python3 /home/arnav.gautam/model/altered_run_model.py $region`)
end

# Collect names and indices of the variables of the economic optimization, which are not named when imported to JuMP
function map_variables_from_AMPL()
    var_name_filename = "/home/arnav.gautam/PowerModelsDistribution.jl/EconomicCostOptimization.col"
    var_name_file = open(var_name_filename, "r")
    Y_a_var_name_to_index = Dict()
    all_var_names_to_index = Dict()

    tech_names = ["C_SOFC", "P_SOFC", "GENERATOR", "PV"]
    tech_dispatch_timeseries_var_indices = Dict(tech_name=>Dict() for tech_name in tech_names)
    let line_number = 1  
        # read till end of file
        while ! eof(var_name_file)
            # read a new / next line for every iteration          
            var_name = readline(var_name_file)
            all_var_names_to_index[var_name] = line_number
            if startswith(var_name, "Y_a[")
                println(var_name)
                Y_a_var_name_to_index[var_name[6:end-2]] = line_number
            elseif startswith(var_name, "X_p[")
                tech_name, timestamp = split(var_name[6:end-1], "',")
                tech_dispatch_timeseries_var_indices[tech_name][timestamp] = line_number
            end
            line_number += 1
        end
    end
    close(var_name_file)
    return (Y_a_var_name_to_index, all_var_names_to_index, tech_dispatch_timeseries_var_indices)
end

# Instantiate the JuMP model from the AMPL economic optimization formulation
function read_economic_model_from_AMPL_into_jump(;my_optimizer = nothing, destination_model = nothing)
    filename = "/home/arnav.gautam/PowerModelsDistribution.jl/EconomicCostOptimization.mps"#.nl"
    if isnothing(destination_model)
        destination_model = JuMP.read_from_file(filename)
    else
        MOI.read_from_file(destination_model, filename)
    end
    if !isnothing(my_optimizer)
        set_optimizer(destination_model, my_optimizer)
    end
    return destination_model
end

# Get the variables associated with the technology numbers from the solved model
function get_tech_number_variables(jump_model, Y_a_var_name_to_index)
    Y_a_variables = Dict()
    all_variables_of_jump_model = JuMP.all_variables(jump_model)
    tech_name_to_var_name = Dict()
    for (tech_name, idx) in Y_a_var_name_to_index
        tech_name_to_var_name[tech_name] = "C$idx"
        Y_a_variables[tech_name] = [var for var in all_variables_of_jump_model if name(var) == tech_name_to_var_name[tech_name]][1]
    end
    # for (tech_name, idx) in Y_a_var_name_to_index
    #     Y_a_variables[tech_name] = VariableRef(jump_model, MOI.VariableIndex(idx))# all_variables_of_jump_model[idx]
    # end
    return Y_a_variables
end

# Extract the technology numbers from the solved model
function get_optimal_system_details(jump_model)
    # Collect variables from AMPL
    Y_a_var_name_to_index, all_var_names_to_index, tech_dispatch_timeseries_var_indices = map_variables_from_AMPL()
    
    # Store references to relevant variables
    Y_a_variables = get_tech_number_variables(jump_model, Y_a_var_name_to_index)
    
    # Compile values for relevant variables
    system_numbers = Dict()
    for tech_name in keys(Y_a_variables)
        system_numbers[tech_name] = JuMP.value(Y_a_variables[tech_name])
    end

    return system_numbers
end

# Make the economic JuMP model into a multi-objective problem (with environmental objectives)
function make_economic_JuMP_model_multi_objective(jump_model, my_inner_optimizer)
    set_optimizer(jump_model, () -> MOA.Optimizer(my_inner_optimizer))
    set_attribute(jump_model, MOA.Algorithm(), MOA.TambyVanderpooten())
    set_attribute(jump_model, MOI.TimeLimitSec(), 36000)

    # Collect variables from AMPL
    Y_a_var_name_to_index, all_var_names_to_index, tech_dispatch_timeseries_var_indices = map_variables_from_AMPL()
    
    # Store references to relevant variables
    Y_a_variables = get_tech_number_variables(jump_model, Y_a_var_name_to_index)

    tech_dispatch_timeseries_variables = Dict(tech_name=>Dict() for tech_name in tech_names)
    for tech_name in keys(tech_dispatch_timeseries_var_indices)
        for (timestamp, idx) in tech_dispatch_timeseries_var_indices[tech_name]
            tech_dispatch_timeseries_variables[tech_name][timestamp] = all_variables_of_jump_model[idx]
        end
    end

    tech_to_CO2 = Dict("C_SOFC"=>0.771618, "P_SOFC"=>0.771618, "GENERATOR"=>1.555, "PV"=>0)

    tech_to_PM = Dict("C_SOFC"=>0, "P_SOFC"=>0, "GENERATOR"=>0.00007, "PV"=>0)

    # Add objective for environment (CO2 emissions)
    CO2_terms = MOI.VectorAffineTerm{Float64}[]
    CO2_output_index = 1
    for tech_name in tech_names
        if tech_to_CO2[tech_name] == 0
            continue
        end
        for timestamp in keys(tech_dispatch_timeseries_variables[tech_name])
            push!(CO2_terms, MOI.VectorAffineTerm(CO2_output_index, MOI.ScalarAffineTerm(tech_to_CO2[tech_name], index(tech_dispatch_timeseries_variables[tech_name][timestamp]))))
        end
    end

    CO2_term_zeros = zeros(length(CO2_terms))

    # Add objective for environment (PM emissions)
    PM_terms = MOI.VectorAffineTerm{Float64}[]
    PM_output_index = 2
    for tech_name in tech_names
        if tech_to_PM[tech_name] == 0
            continue
        end
        for timestamp in keys(tech_dispatch_timeseries_variables[tech_name])
            push!(PM_terms, MOI.VectorAffineTerm(PM_output_index, MOI.ScalarAffineTerm(tech_to_PM[tech_name], index(tech_dispatch_timeseries_variables[tech_name][timestamp]))))
        end
    end

    PM_term_zeros = zeros(length(PM_terms))

    economic_terms = MOI.VectorAffineTerm{Float64}[]
    economic_output_index = 3
    for (econ_var, econ_coeff) in pairs(economic_objective.terms)
        push!(economic_terms, MOI.VectorAffineTerm(economic_output_index, MOI.ScalarAffineTerm(econ_coeff, index(econ_var))))
    end
    economic_term_zeros = zeros(length(economic_terms))

    all_terms = [CO2_terms ; PM_terms ; economic_terms]
    all_zeros = [CO2_term_zeros ; PM_term_zeros ; economic_term_zeros]

    multi_objective_function = MOI.VectorAffineFunction(
            all_terms,
            all_zeros
        )

    MOI.set(jump_model, MOI.ObjectiveFunction{typeof(multi_objective_function)}(), multi_objective_function)
    return jump_model
end

# Create an BilevelJuMP model containing unlinked economic optimization and TPIA
function construct_bilevel_jump_model(my_optimizer)
    bi_level_jump_model = BilevelModel(my_optimizer)
    set_attribute(bi_level_jump_model.upper, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    set_attribute(bi_level_jump_model.lower, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    read_economic_model_from_AMPL_into_jump(;destination_model=backend(bi_level_jump_model.upper));
    create_PMD_JuMP_model(sf_math, LinDist3FlowPowerModel, build_mn_mc_tpia_L1, my_baron_optimizer; PMD_jump_model = backend(bi_level_jump_model.lower));
    # set_optimizer_attribute(bi_level_jump_model, "NonConvex", 2) # If I have Nonconvex/nonconcave/nonpsd objective/constraint error in a MIP solver like Gurobi.
    return bi_level_jump_model
end

# Make the upper level of a BilevelJump model into a multi-objective optimization
function make_bilevel_JuMP_model_multi_objective(bi_level_jump_model, my_inner_optimizer)
    set_optimizer(bi_level_jump_model, () -> MOA.Optimizer(my_inner_optimizer))
    set_attribute(bi_level_jump_model, MOA.Algorithm(), MOA.TambyVanderpooten())
    set_attribute(bi_level_jump_model, MOI.TimeLimitSec(), 36000)
    make_economic_JuMP_model_multi_objective(backend(bi_level_jump_model.upper), my_inner_optimizer)
end

# Change the multi-objective upper problem to use the lower variables of power delivered for tech selection/evaluation
function connect_bi_level_model(bi_level_jump_model)
    # Not implemented yet
end

# Change the lower problem to incorporate the upper variables
function make_bilevel_JuMP_model_equity_aware(bi_level_jump_model)
    # Not implemented yet
end

#################################################################################################
# HELPER FUNCTIONS                                                                              #
#################################################################################################

function create_sum_of_loads_DAT_file()
    # Locate the loads of the region
    loads_filename = joinpath(FILENAME_BASE, "Loads.dss")
    loadshape_base_values = Dict()
    loads_file = open(loads_filename, "r")
    while !eof(loads_file)
        load_line = readline(loads_file)
        if length(load_line) == 0
            continue
        end
        load_array = split(load_line, " ")
        load_p_base_value = parse(Float64, load_array[9][4:end])
        load_q_base_value = parse(Float64, load_array[10][6:end])
        loadshape_name = load_array[12][8:end]
        loadshape_base_values[loadshape_name] = (load_p_base_value, load_q_base_value)
    end
    close(loads_file)


    # Locate the loadshapes of the region
    loadshape_reference_filename = joinpath(FILENAME_BASE, "LoadShapes.dss")
    loadshape_file_list = Dict()
    loadshape_reference_file = open(loadshape_reference_filename, "r")
    while !eof(loadshape_reference_file)
        loadshape_reference_line = readline(loadshape_reference_file)
        if length(loadshape_reference_line) == 0
            continue
        end
        loadshape_reference_array = split(loadshape_reference_line, " ")
        loadshape_name = loadshape_reference_array[2][11:end]
        mult_path = loadshape_reference_array[8][7:end-1]
        qmult_path = loadshape_reference_array[11][7:end-1]
        loadshape_file_list[loadshape_name] = (mult_path, qmult_path)
    end
    close(loadshape_reference_file)

    # Read in the files to a format Julia understands
    # For each timestamp, calculate total load
    power_timeseries = nothing
    for (loadshape_name, (p_timeseries_file, q_timeseries_file)) in pairs(loadshape_file_list)
        load_p_base_value, load_q_base_value = loadshape_base_values[loadshape_name]
        p_timeseries = readdlm(normpath(joinpath(loadshape_reference_filename, "..", p_timeseries_file)))
        q_timeseries = readdlm(normpath(joinpath(loadshape_reference_filename, "..", q_timeseries_file)))
        new_power_timeseries = ((load_p_base_value.*p_timeseries).^2 + (load_q_base_value.*q_timeseries).^2).^0.5
        if isnothing(power_timeseries)
            power_timeseries = new_power_timeseries
        else
            power_timeseries += new_power_timeseries
        end
    end

    # Write out DAT file with data
    return write_DAT_file(power_timeseries; DAT_file_name_suffix="_sum_of")
end

function ensure_tpia_output_is_nonzero(slack_power_filename)
    total_slack_power = 0
    slack_power_file = open(slack_power_filename, "r")
    for slack_power_line in readlines(slack_power_file)
        line_values = split(slack_power_line, " ", keepempty=false)
        if length(line_values) == 3 && line_values[2] == "total"
            total_slack_power += parse(Float64, line_values[3])
        end
    end
    close(slack_power_file)
    @assert(total_slack_power != 0, "Slack power is zero")
end

function run_economic_opti_on_text(; DAT_file_name_suffix="")
    refresh_economic_model_in_AMPL("$SFO_SMART_DS_REGION$DAT_file_name_suffix")
    economic_jump_model = read_economic_model_from_AMPL_into_jump(;my_optimizer=my_gurobi_optimizer)
    JuMP.optimize!(economic_jump_model)
    return economic_jump_model
end

# Give a full summary of a solved bilevel optimization
function summarize_bilevel_opti_solution(bi_level_jump_model)
    JuMP.solution_summary(bi_level_jump_model)
    primal_status(bi_level_jump_model)
    objective_value(bi_level_jump_model)
    lower_obj_type = JuMP.objective_function_type(Lower(bi_level_jump_model))
    JuMP.objective_function(Lower(bi_level_jump_model), lower_obj_type)
    upper_obj_type = JuMP.objective_function_type(Upper(bi_level_jump_model))
    JuMP.objective_function(Upper(bi_level_jump_model), upper_obj_type)
end

function make_all_timeseries_files_multi_year()
    # Read in single_year_timeseries_files
    single_year_timeseries_files = []

    # Concat each one to create a longer timeseries

    # Set the new timeseries as the contents of the files

    # Return the old file contents to be restored later
    return single_year_timeseries_files
end

function reset_all_timeseries_files(single_year_timeseries_files)
    # Restore single_year_timeseries_files
end

#################################################################################################
# CONFIGURATION 1: Economic Opti(sum(loads))                                                    #
#################################################################################################
# create_sum_of_loads_DAT_file()
# economic_jump_model = run_economic_opti_on_text(; DAT_file_name_suffix="_sum_of")
# optimal_system_details = get_optimal_system_details(economic_jump_model)
# obj_value = objective_value(economic_jump_model)
# JuMP.solution_summary(economic_jump_model)

# TODO add functionality to copy relevant files to be saved for later analysis

#################################################################################################
# CONFIGURATION 2.1: TMPIA(single timestamp) → text output → Economic Opti(text output)         #
#################################################################################################
sf_math_timestamp = create_PMD_mathematical_model(sf_casefile; time_series_to_run="daily")
sf_math_timestamp["nw"]["0"]["gen"]["1"]["pmax"] = [0.0,0.0,0.0]
sf_math_timestamp["nw"]["0"]["gen"]["1"]["pmin"] = [0.0,0.0,0.0]
sf_math_timestamp["nw"]["0"]["gen"]["1"]["qmax"] = [0.0,0.0,0.0]
sf_math_timestamp["nw"]["0"]["gen"]["1"]["qmin"] = [0.0,0.0,0.0]
pmd_timestamp_ACPU_jump_model = create_PMD_JuMP_model(sf_math_timestamp, ACPUPowerModel, build_mc_tmpia, tighter_ipopt_optimizer)
JuMP.optimize!(pmd_timestamp_ACPU_jump_model)
total_pg = [JuMP.value(var) for var in JuMP.all_variables(pmd_timestamp_ACPU_jump_model) if occursin("_pg_", name(var))]
total_qg = [JuMP.value(var) for var in JuMP.all_variables(pmd_timestamp_ACPU_jump_model) if occursin("_qg_", name(var))]
total_power_generated = (total_pg.^2 + total_qg.^2).^0.5
p_slack_cap_dict = Dict([(name(var),JuMP.value(var)) for var in JuMP.all_variables(pmd_timestamp_ACPU_jump_model) if occursin("_p_slack_cap_", name(var))])
q_slack_cap_dict = Dict([(name(var),JuMP.value(var)) for var in JuMP.all_variables(pmd_timestamp_ACPU_jump_model) if occursin("_q_slack_cap_", name(var))])
# total_slack_power = (p_slack_cap.^2 + q_slack_cap.^2).^0.5

# slack_power_output_map = Dict()
slack_power_df = DataFrame(bus_idx = Any[], slack_power = Float64[], phase_str = String[])
for bus_idx in keys(sf_math_timestamp["nw"]["0"]["bus"])
    p_slack_cap_val = 0
    q_slack_cap_val = 0
    multi_phase = true
    phase_str = ""
    for phase in ["[1]", "[2]", "[3]"]
        p_slack_cap_var_name = string("0_p_slack_cap_", bus_idx, phase)
        if haskey(p_slack_cap_dict, p_slack_cap_var_name)
            phase_str *= phase
            p_slack_cap_val += p_slack_cap_dict[p_slack_cap_var_name]
        else
            multi_phase = false
        end
        q_slack_cap_var_name = string("0_q_slack_cap_", bus_idx, phase)
        if haskey(q_slack_cap_dict, q_slack_cap_var_name)
            q_slack_cap_val += q_slack_cap_dict[q_slack_cap_var_name]
        else
            multi_phase = false
        end
    end
    total_slack_power = (p_slack_cap_val.^2 + q_slack_cap_val.^2).^0.5
    # slack_power_output_map[bus_idx] = total_slack_power
    push!(slack_power_df, (bus_idx = bus_idx, slack_power = total_slack_power, phase_str = phase_str))
end

CSV.write("slack_power_output.csv", slack_power_df)

# #################################################################################################
# # CONFIGURATION 2: TPIA(single timestamp) → text output → Economic Opti(text output)            #
# #################################################################################################
# sf_math_timestamp = create_PMD_mathematical_model(sf_casefile; time_series_to_run="daily")

# # Solve TPIA, either:

#     # DIRECTLY
#         # tpia_model = create_PMD_JuMP_model(sf_math_timestamp, LinDist3FlowPowerModel, build_mn_mc_tpia_L1, tighter_ipopt_optimizer)
#         # JuMP.optimize!(tpia_model)
#         # JuMP.solution_summary(tpia_model)

#     # OR

#     # DEBUGGING
#         # pmd_timestamp_LinDist3Flow_jump_model = create_PMD_JuMP_model(sf_math_timestamp, LinDist3FlowPowerModel, build_mn_mc_tpia_L1, my_baron_optimizer)
#         # JuMP.optimize!(pmd_timestamp_LinDist3Flow_jump_model)
#         # total_pg = [JuMP.value(var) for var in JuMP.all_variables(pmd_timestamp_LinDist3Flow_jump_model) if occursin("_pg_", name(var))]
#         # total_qg = [JuMP.value(var) for var in JuMP.all_variables(pmd_timestamp_LinDist3Flow_jump_model) if occursin("_qg_", name(var))]
#         # total_power_generated = (total_pg.^2 + total_qg.^2).^0.5

#         pmd_timestamp_ACPU_jump_model = create_PMD_JuMP_model(sf_math_timestamp, ACPUPowerModel, build_mn_mc_tpia_L1, tighter_ipopt_optimizer)
#         JuMP.optimize!(pmd_timestamp_ACPU_jump_model)
#         total_pg = [JuMP.value(var) for var in JuMP.all_variables(pmd_timestamp_ACPU_jump_model) if occursin("_pg_", name(var))]
#         total_qg = [JuMP.value(var) for var in JuMP.all_variables(pmd_timestamp_ACPU_jump_model) if occursin("_qg_", name(var))]
#         total_power_generated = (total_pg.^2 + total_qg.^2).^0.5
#         total_p_slack = [JuMP.value(var) for var in JuMP.all_variables(pmd_timestamp_ACPU_jump_model) if occursin("_p_slack_in_", name(var))]

# LinDist3Flow_slack_power_timeseries_data = run_PMD_directly_and_write_result_to_file(sf_math_timestamp, LinDist3FlowPowerModel, solve_mn_mc_tpia_L1, my_baron_optimizer; is_multinetwork=true, write_to_DAT_file=false)
# @assert(sum(LinDist3Flow_slack_power_timeseries_data) != 0, "Slack power is zero")

# ACPU_slack_power_timeseries_data = run_PMD_directly_and_write_result_to_file(sf_math_timestamp, ACPUPowerModel, solve_mn_mc_tpia_L1, my_baron_optimizer; is_multinetwork=true, write_to_DAT_file=false)
# @assert(sum(ACPU_slack_power_timeseries_data) != 0, "Slack power is zero")


# # ensure_tpia_output_is_nonzero(pmd_timestamp_result_filename)

# economic_jump_model = run_economic_opti_on_text()
# optimal_system_details = get_optimal_system_details(economic_jump_model)
# obj_value = objective_value(economic_jump_model)
# # JuMP.solution_summary(economic_jump_model)

# # TODO add functionality to copy relevant files to be saved for later analysis

# #################################################################################################
# # CONFIGURATION 3: TPIA(one-year timestamp) → text output → Economic Opti(text output)          #
# #################################################################################################
# sf_math = create_PMD_mathematical_model(sf_casefile; time_series_to_run="yearly")
# pmd_timeseries_result_filename = run_PMD_directly_and_write_result_to_file(sf_math, LinDist3FlowPowerModel, solve_mn_mc_tpia_L1, my_baron_optimizer; is_multinetwork=true)
# ensure_tpia_output_is_nonzero(pmd_timeseries_result_filename)

# economic_jump_model = run_economic_opti_on_text()
# # JuMP.solution_summary(economic_jump_model)

# # TODO add functionality to copy relevant files to be saved for later analysis

# #################################################################################################
# # CONFIGURATION 4: TPIA(one-year timestamp) → Economic Opti(TPIA result variables)              #
# #################################################################################################
# bi_level_jump_model = construct_bilevel_jump_model(my_baron_optimizer)
# JuMP.optimize!(bi_level_jump_model)
# # summarize_bilevel_opti_solution(bi_level_jump_model)

# connect_bi_level_model(bi_level_jump_model)
# JuMP.optimize!(bi_level_jump_model)
# # summarize_bilevel_opti_solution(bi_level_jump_model)

# # TODO add functionality to copy relevant files to be saved for later analysis

# #################################################################################################
# # CONFIGURATION 5: Equitable TPIA(one-year timestamp) → Economic Opti(TPIA result variables)    #
# #################################################################################################
# make_bilevel_JuMP_model_equity_aware(bi_level_jump_model)
# JuMP.optimize!(bi_level_jump_model)
# # JuMP.solution_summary(MOP_jump_model)

# # TODO add functionality to copy relevant files to be saved for later analysis

# #################################################################################################
# # CONFIGURATION 6: TPIA(one-year timestamp) → Multi Objective Opti(TPIA result variables)       #
# #################################################################################################
# bi_level_jump_model = construct_bilevel_jump_model(my_baron_optimizer)
# make_bilevel_JuMP_model_multi_objective(bi_level_jump_model, my_baron_optimizer)
# JuMP.optimize!(bi_level_jump_model)
# # JuMP.solution_summary(MOP_jump_model)

# # TODO add functionality to copy relevant files to be saved for later analysis

# #################################################################################################
# # CONFIGURATION 7: TPIA(multi-year timestamp) → Multi Objective Opti(TPIA result variables)     #
# #################################################################################################
# single_year_timeseries_files = make_all_timeseries_files_multi_year()
# sf_math = create_PMD_mathematical_model(sf_casefile; time_series_to_run="yearly")
# bi_level_jump_model = construct_bilevel_jump_model(my_baron_optimizer)
# make_bilevel_JuMP_model_multi_objective(bi_level_jump_model, my_baron_optimizer)
# JuMP.optimize!(bi_level_jump_model)
# # JuMP.solution_summary(bi_level_jump_model)
# reset_all_timeseries_files(single_year_timeseries_files)

# # TODO add functionality to copy relevant files to be saved for later analysis