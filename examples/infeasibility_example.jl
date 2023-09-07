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
# using PyCall
# using AmplNLWriter

# function Couenne_Optimizer()
#     return AmplNLWriter.Optimizer("couenne")
# end

# Two ways to obtain convergence
# my_ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>200, "tol"=>100.0, "dual_inf_tol"=>100.0, "constr_viol_tol"=>1.0, "compl_inf_tol"=>0.01)
my_ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>5000, "acceptable_tol"=>100000.0, "acceptable_dual_inf_tol"=>10000.0, "acceptable_constr_viol_tol"=>1.0, "acceptable_compl_inf_tol"=>1.0, "acceptable_obj_change_tol"=>100.0)
tighter_ipopt_optimizer = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>5000, "acceptable_tol"=>100000.0, "acceptable_dual_inf_tol"=>10000.0, "acceptable_constr_viol_tol"=>1.0, "acceptable_compl_inf_tol"=>1.0, "acceptable_obj_change_tol"=>10.0)
my_gurobi_optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit"=>36000, "MipGap"=>0.05, "nonconvex"=>2) # "gurobi_options"=>"timelim 300 mipgap 0.05 nonconvex 2"
my_baron_optimizer = optimizer_with_attributes(
    BARON.Optimizer,
    "threads"=>16,
    # "EpsA"=>1e6,
    # "EpsR"=>1e6,
    # "DeltaTerm"=>1,
    # "DeltaT"=>600,
    # "AbsConFeasTol"=>1e6,
    # "RelConFeasTol"=>0.1,
    # "AbsIntFeasTol"=>1e6,
    # "RelIntFeasTol"=>0.1,
    # "FirstFeas"=>1,
    # "FirstLoc"=>1,
    "MaxTime"=>36000
    )

sfo_smart_ds_region = "P4U"

# sf_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_$(sfo_smart_ds_region)_Timeseries/scenarios/base_timeseries/opendss/Master.dss"
sf_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_$(sfo_smart_ds_region)_Timeseries/scenarios/base_timeseries/opendss/p4uhs0_4/p4uhs0_4--p4udt0/Master.dss"
sf_eng = parse_file(sf_casefile; multinetwork=true, time_series="yearly")
set_time_elapsed!(sf_eng, 0.25)

# TODO: find a way to remove the infinite power of the slack bus
# sf_eng["voltage_source"]["source"]["status"]=DISABLED
# sf_eng["voltage_source"]["source"]["rs"] .= 1e100#Inf
# sf_eng["voltage_source"]["source"]["xs"] .= 1e100#Inf
# sf_eng["voltage_source"]["source"]["pg_ub"] = [0, 0, 0]
# sf_eng["voltage_source"]["source"]["qg_ub"] = [0, 0, 0]

# success?
# for nw in 1:length(sf_eng["nw"])
#     sf_eng["nw"]["$(nw-1)"]["voltage_source"]["source"]["pg_ub"] = [0, 0, 0]
#     sf_eng["nw"]["$(nw-1)"]["voltage_source"]["source"]["qg_ub"] = [0, 0, 0]
# end

passthrough_dict =  Dict{String,Vector{String}}(
    "bus" => String["lat", "lon"])
sf_math = transform_data_model(sf_eng, eng2math_passthrough=passthrough_dict)

# for nw in 1:length(sf_eng["nw"])
#     @assert(sf_math["nw"]["$(nw-1)"]["gen"]["1"]["pmax"] == [0.0,0.0,0.0])
#     @assert(sf_math["nw"]["$(nw-1)"]["gen"]["1"]["qmax"] == [0.0,0.0,0.0])
# end

tpia_model = JuMP.Model()
instantiate_mc_model(sf_math, LinDist3FlowPowerModel, build_mn_mc_tpia_L1; jump_model = tpia_model);
set_optimizer(tpia_model, tighter_ipopt_optimizer)
JuMP.optimize!(tpia_model)
JuMP.solution_summary(tpia_model)


# if false

sf_tpia_solution = solve_mn_mc_tpia_L1(sf_math, LinDist3FlowPowerModel, tighter_ipopt_optimizer; multinetwork=true)

slack_power_file = open("/home/arnav.gautam/model/data/$(sfo_smart_ds_region)_load.dat", "w")
write(slack_power_file, "param: d_p1 := \n")
for nw_num in 1:length(keys(sf_tpia_solution["solution"]["nw"]))
    solution_p_slack = sum(sum(val["p_slack_out"]) - sum(val["p_slack_in"]) for val in values(sf_tpia_solution["solution"]["nw"]["$(nw_num-1)"]["bus"]))
    solution_q_slack = sum(sum(val["q_slack_out"]) - sum(val["q_slack_in"]) for val in values(sf_tpia_solution["solution"]["nw"]["$(nw_num-1)"]["bus"]))
    solution_slack_power = sqrt(solution_p_slack^2 + solution_q_slack^2)
    write(slack_power_file, "$nw_num  total $solution_slack_power\n")
    write(slack_power_file, "$nw_num  cl $solution_slack_power\n")
    write(slack_power_file, "$nw_num  rl 0.0\n")
end
write(slack_power_file, ";\n")
close(slack_power_file)

# end

# if false

#     richmond_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_P19U/Master.dss"
#     # sf_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_P4U/Master.dss"
#     # sf_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_P4U_Timeseries/scenarios/base_timeseries/opendss/p4uhs0_4/Master.dss"
#     richmond_eng = parse_file(richmond_casefile)

#     richmond_eng["voltage_source"]["source"]["rs"] .= Float64.Inf
#     richmond_eng["voltage_source"]["source"]["xs"] .= Float64.Inf


#     richmond_math = transform_data_model(richmond_eng, eng2math_passthrough=passthrough_dict)

#     richmond_tpia_solution = solve_mc_tpia_L1(richmond_math, ACPUPowerModel, my_ipopt_optimizer)
#     solution_p_slack = sum(sum(val["p_slack_in"]) - sum(val["p_slack_out"]) for val in values(richmond_tpia_solution["solution"]["bus"]))
#     solution_q_slack = sum(sum(val["q_slack_in"]) - sum(val["q_slack_out"]) for val in values(richmond_tpia_solution["solution"]["bus"]))
#     solution_slack_power = sqrt(solution_p_slack^2 + solution_q_slack^2)

# end

# py"""
# import os
# run_model_path = os.path.join(os.path.dirname(os.getcwd()), "model", "run_model.py")
# with open(run_model_path) as f:
#     exec(run_model_path.read())
# """

# if false

run(`python3 /home/arnav.gautam/model/altered_run_model.py $sfo_smart_ds_region`)

# end
# run(`python3 /home/arnav.gautam/model/altered_run_model.py R1-1247-2`)

# if false
# INSTANTIATE THE JuMP MODEL FROM JAIME'S OPTIMIZATION FORMULATION
filename = "/home/arnav.gautam/PowerModelsDistribution.jl/EconomicCostOptimization.mps"#.nl"
var_name_filename = "/home/arnav.gautam/PowerModelsDistribution.jl/EconomicCostOptimization.col"
var_name_file = open(var_name_filename, "r")

tech_names = ["C_SOFC", "P_SOFC", "GENERATOR", "PV"]

Y_a_var_name_to_index = Dict()
all_var_names_to_index = Dict()
tech_dispatch_timeseries_var_indices = Dict(tech_name=>Dict() for tech_name in tech_names)
let line_number = 1  
    # read till end of file
    while ! eof(var_name_file)
        # read a new / next line for every iteration          
        var_name = readline(var_name_file)
        all_var_names_to_index[var_name] = line_number
        if startswith(var_name, "Y_a[")
            Y_a_var_name_to_index[var_name[6:end-2]] = line_number
        elseif startswith(var_name, "X_p[")
            tech_name, timestamp = split(var_name[6:end-1], "',")
            tech_dispatch_timeseries_var_indices[tech_name][timestamp] = line_number
        end
        line_number += 1
    end
end
close(var_name_file)

# constraint_obj_name_filename = "/home/arnav.gautam/model/EconomicCostOptimization.row"
# src =
#     MOI.FileFormats.Model(; format = MOI.FileFormats.FORMAT_NL, filename = filename, kwargs...)
# MOI.read_from_file(src, filename)
# MOI.copy_to(Upper(bi_level_jump_model), src)
economic_jump_model = JuMP.read_from_file(filename)

economic_objective = objective_function(economic_jump_model);

set_optimizer(economic_jump_model, my_gurobi_optimizer)

JuMP.optimize!(economic_jump_model)

JuMP.solution_summary(economic_jump_model)

MOP_jump_model = JuMP.read_from_file(filename)

# Make this a multi-objective problem

set_optimizer(MOP_jump_model, () -> MOA.Optimizer(my_gurobi_optimizer)) # my_baron_optimizer) # my_ipopt_optimizer)#
set_attribute(MOP_jump_model, MOA.Algorithm(), MOA.TambyVanderpooten())
# set_attribute(MOP_jump_model, MOA.SolutionLimit(), 10)
set_attribute(MOP_jump_model, MOI.TimeLimitSec(), 36000)

# Store references to relevant variables
Y_a_variables = Dict()
all_variables_of_MOP_jump_model = JuMP.all_variables(MOP_jump_model)
for (tech_name, idx) in Y_a_var_name_to_index
    Y_a_variables[tech_name] = all_variables_of_MOP_jump_model[idx]
end

tech_dispatch_timeseries_variables = Dict(tech_name=>Dict() for tech_name in tech_names)
for tech_name in keys(tech_dispatch_timeseries_var_indices)
    for (timestamp, idx) in tech_dispatch_timeseries_var_indices[tech_name]
        tech_dispatch_timeseries_variables[tech_name][timestamp] = all_variables_of_MOP_jump_model[idx]
    end
end


tech_to_CO2 = Dict("C_SOFC"=>0.771618, "P_SOFC"=>0.771618, "GENERATOR"=>1.555, "PV"=>0)

tech_to_PM = Dict("C_SOFC"=>0, "P_SOFC"=>0, "GENERATOR"=>0.00007, "PV"=>0)

# # Add environment impact parameters (CO2 emissions)
# CO2_variables = Dict()
# # co_constraints = Dict()
# # all_variables_of_MOP_jump_model = JuMP.all_variables(MOP_jump_model)

# # Add constant for CO2 environmental impact (JK not doing it this way because it makes the objective quadratic)
# for (tech_name, var) in Y_a_variables
#     # println(tech_name)
#     # tech_num_var = all_variables_of_MOP_jump_model[idx]
#     # Create new variable for CO2 environmental impact of this tech
#     tech_CO2_var = JuMP.@variable(MOP_jump_model, base_name="$(tech_name)_CO2")
#     CO2_variables[tech_name] = tech_CO2_var
#     # Fix variable at known amount
#     fix(tech_CO2_var, tech_to_CO2[tech_name])
# end

# Add objective for environment (CO2 emissions)


CO2_terms = MOI.VectorAffineTerm{Float64}[]
CO2_output_index = 1
# let i = 1
for tech_name in tech_names
    if tech_to_CO2[tech_name] == 0
        continue
    end
    for timestamp in keys(tech_dispatch_timeseries_variables[tech_name])
        push!(CO2_terms, MOI.VectorAffineTerm(CO2_output_index, MOI.ScalarAffineTerm(tech_to_CO2[tech_name], index(tech_dispatch_timeseries_variables[tech_name][timestamp]))))
    end
end
# end
CO2_term_zeros = zeros(length(CO2_terms))
# CO2_objective = MOI.VectorAffineFunction(CO2_terms, CO2_term_zeros)#JuMP.@expression(MOP_jump_model, # Min,


# CO2_objective = JuMP.@expression(MOP_jump_model, # Min,
#     sum(
#         sum(
#             tech_dispatch_timeseries_variables[tech_name][timestamp] for timestamp in keys(tech_dispatch_timeseries_variables[tech_name])
#         ) * tech_to_CO2[tech_name] for tech_name in tech_names) # * Y_a_variables[tech_name] (this is number procured, not number running. and I think it's accounted for in X_p)
# )

# # Add environment impact parameters (PM emissions)
# PM_variables = Dict()
# # co_constraints = Dict()
# # all_variables_of_MOP_jump_model = JuMP.all_variables(MOP_jump_model)
# for (var_name, var) in Y_a_variables
#     # println(var_name)
#     # tech_num_var = all_variables_of_MOP_jump_model[idx]
#     # Create new variable for PM environmental impact of this tech
#     tech_PM_var = JuMP.@variable(MOP_jump_model, base_name="$(var_name)_PM")
#     PM_variables[var_name] = tech_PM_var
#     # Fix variable at known amount
#     fix(tech_PM_var, tech_to_PM[var_name])
# end

# for tech_name in tech_names
#     println(tech_to_PM[tech_name])
#     # for timestamp in keys(tech_dispatch_timeseries_variables[tech_name])
#     #     println(index(tech_dispatch_timeseries_variables[tech_name][timestamp]))
#     #     println(tech_dispatch_timeseries_variables[tech_name][timestamp])
#     # end
# end


# Add objective for environment (PM emissions)
PM_terms = MOI.VectorAffineTerm{Float64}[]
PM_output_index = 2
# let i = 1
for tech_name in tech_names
    if tech_to_PM[tech_name] == 0
        continue
    end
    for timestamp in keys(tech_dispatch_timeseries_variables[tech_name])
        push!(PM_terms, MOI.VectorAffineTerm(PM_output_index, MOI.ScalarAffineTerm(tech_to_PM[tech_name], index(tech_dispatch_timeseries_variables[tech_name][timestamp]))))
    end
end
# end
PM_term_zeros = zeros(length(PM_terms))
# PM_objective = MOI.VectorAffineFunction(PM_terms, PM_term_zeros)#JuMP.@expression(MOP_jump_model, # Min,
#     [
#         (
#             MOI.VectorAffineTerm(i, MOI.ScalarAffineTerm(tech_to_PM[tech_name], tech_dispatch_timeseries_variables[tech_name][timestamp])),
#             1
#         ) for timestamp in keys(tech_dispatch_timeseries_variables[tech_name]) for tech_name in tech_names
#     ]
# )

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
        # fill(0.0, 3),
    )


MOI.set(MOP_jump_model, MOI.ObjectiveFunction{typeof(multi_objective_function)}(), multi_objective_function)

JuMP.optimize!(MOP_jump_model)

# TODO RUN MOP_jump_model
# end
# economic_variables = all_variables(MOP_jump_model)
# economic_constraints = Vector{Any}()
# println(list_of_constraint_types(MOP_jump_model))
# for (function_type, set_type) in list_of_constraint_types(MOP_jump_model)
#     println(all_constraints(MOP_jump_model, function_type, set_type))
# end

JuMP.solution_summary(MOP_jump_model)

# for (var, idx) in pairs(var_name_to_index)
#     # println(JuMP.variable_by_name(MOP_jump_model, "C$idx"))
#     # println(JuMP.value(JuMP.variable_by_name(MOP_jump_model, "C$idx")))
#     println(JuMP.all_variables(MOP_jump_model)[idx])
#     println(JuMP.value(JuMP.all_variables(MOP_jump_model)[idx]))
#     # println("$(var) = $(JuMP.value(JuMP.variable_by_name(MOP_jump_model, "C$idx")))")#JuMP.all_variables(MOP_jump_model)[idx]))")
# end

# all_variables_of_MOP_jump_model = JuMP.all_variables(MOP_jump_model)
# for (var_name, idx) in Y_a_var_name_to_index
#     println(var_name)
#     var = all_variables_of_MOP_jump_model[idx]
#     println(JuMP.value(var))
# end

# for (idx, var) in pairs(JuMP.all_variables(MOP_jump_model))
#     if (all_indices_to_var_name[idx] in keys(var_name_to_index))
#         println(var)
#         println(all_vars_to_index[idx])
#         println(JuMP.value(var))
#     end
# end

# JuMP.value(JuMP.all_variables(MOP_jump_model))

bi_level_jump_model = BilevelModel(() -> MOA.Optimizer(my_baron_optimizer))
set_attribute(bi_level_jump_model, MOA.Algorithm(), MOA.TambyVanderpooten())
set_attribute(bi_level_jump_model, MOI.TimeLimitSec(), 36000)

# JuMP.set_optimizer(bi_level_jump_model.upper, () -> MOA.Optimizer(my_baron_optimizer))


set_optimizer(bi_level_jump_model.upper, my_gurobi_optimizer)
MOI.read_from_file(backend(bi_level_jump_model.upper), filename) # get MOP_jump_model into bilevel
set_attribute(bi_level_jump_model.upper, MOI.ObjectiveSense(), MOI.MIN_SENSE)

# MOI.copy_to(backend(bi_level_jump_model.upper), backend(MOP_jump_model))

# Change the multi-objective upper problem to use the lower variables of power delivered for tech selection/evaluation


# test_jump_model = JuMP.Model()

# tpia_im_model = 
set_optimizer(bi_level_jump_model.lower, my_ipopt_optimizer)
instantiate_mc_model(sf_math, LinDist3FlowPowerModel, build_mc_tpia_L1; jump_model = bi_level_jump_model.lower); # test_jump_model) # ACPUPowerModel, OR LinDist3FlowPowerModel, OR FBSUBFPowerModel, OR FOTPUPowerModel, OR FOTRUPowerModel
set_attribute(bi_level_jump_model.lower, MOI.ObjectiveSense(), MOI.MIN_SENSE)

# tpia_jump_model = tpia_im_model.model
# tpia_variables = all_variables(tpia_jump_model)
# tpia_constraints = Vector{Any}()
# for (function_type, set_type) in list_of_constraint_types(tpia_jump_model)
#     append!(tpia_constraints, all_constraints(tpia_jump_model, function_type, set_type))
# end

# set_optimizer_attribute(bi_level_jump_model, "NonConvex", 2) # If I have Nonconvex/nonconcave/nonpsd objective/constraint error in a MIP solver like Gurobi.

# @variable(Lower(bi_level_jump_model), tpia_variables[1:length(tpia_variables)])

# for upper_var in jaime_variables
#     base_name_value = name(upper_var)
#     @variable(Upper(bi_level_jump_model), name(upper_var))
# end

# @variable(Upper(bi_level_jump_model), jaime_variables[1:length(jaime_variables)])

# set_optimizer(bi_level_jump_model, )

# for (function_type, set_type) in tpia_constrai

# Solve the bilevel optimization problem
JuMP.optimize!(bi_level_jump_model)

primal_status(bi_level_jump_model)
objective_value(bi_level_jump_model)
lower_obj_type = JuMP.objective_function_type(Lower(bi_level_jump_model))
JuMP.objective_function(Lower(bi_level_jump_model), lower_obj_type)
upper_obj_type = JuMP.objective_function_type(Upper(bi_level_jump_model))
JuMP.objective_function(Upper(bi_level_jump_model), upper_obj_type)