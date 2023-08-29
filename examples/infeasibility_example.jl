using Pkg
Pkg.activate(pwd())
using Revise
using PowerModelsDistribution
using Ipopt
using BARON
using Gurobi
using JuMP
using BilevelJuMP
# using AmplNLWriter
# using AmplNLReader
using Infiltrator
# import MathOptInterface as MOI
import MultiObjectiveAlgorithms as MOA

# richmond_casefile = "/mnt/c/Users/Arnav/Documents/Research/Data/SMART-DS/SFO_Regions/SFO/P19U/Master.dss"
sf_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/arnav_data/SFO_P4U/Master.dss"

# example_casefile = "/home/arnav.gautam/PowerModelsDistribution.jl/test/data/opendss/case3_balanced.dss"

# richmond_eng = parse_file(richmond_casefile)
sf_eng = parse_file(sf_casefile)

# example_eng = parse_file(example_casefile)

# example_tpia_soln = solve_mc_tpia(example_eng, ACPUPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>500, "tol"=>2))

# function Couenne_Optimizer()
#     return AmplNLWriter.Optimizer("couenne")
# end


# function Ipopt_Optimizer()
#     return AmplNLWriter.Optimizer("ipopt")
# end

# sf_acp_solution = solve_mc_pf(sf_eng, ACPUPowerModel, Ipopt.Optimizer)
# sf_ivr_solution = solve_mc_pf(sf_eng, IVRUPowerModel, Ipopt.Optimizer)
# example_acp_solution = solve_mc_pf(example_eng, ACPUPowerModel, Ipopt.Optimizer)
# example_ivr_solution = solve_mc_pf(example_eng, IVRUPowerModel, Ipopt.Optimizer)


# richmond_solution = solve_mc_pf(richmond_eng, ACPUPowerModel, BARON.Optimizer)
# richmond_solution = solve_mc_pf_pbs(richmond_eng, ACPUPowerModel, Ipopt.Optimizer)

# richmond_solution = solve_mc_tpia(richmond_eng, ACPUPowerModel, Ipopt.Optimizer)
# sf_tpia_solution = solve_mc_tpia(sf_eng, ACPUPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>500, "tol"=>2.0))

# richmond_solution = solve_mc_tpia(richmond_eng, ACPUPowerModel, BARON.Optimizer)

# example_solution = solve_mc_pf(example_eng, ACPUPowerModel, Ipopt.Optimizer)
# example_solution = solve_mc_pf(example_eng, ACPUPowerModel, Couenne_Optimizer)

passthrough_dict =  Dict{String,Vector{String}}(
    "bus" => String["lat", "lon"])
# example_math = transform_data_model(example_eng, eng2math_passthrough=passthrough_dict)
sf_math = transform_data_model(sf_eng, eng2math_passthrough=passthrough_dict)

# ARNAV_solver = optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>200, "acceptable_tol"=>100.0, "acceptable_dual_inf_tol"=>10000.0, "acceptable_constr_viol_tol"=>1.0, "acceptable_compl_inf_tol"=>0.1, "acceptable_obj_change_tol"=>10.0)

bi_level_jump_model = BilevelModel(() -> MOA.Optimizer(Gurobi.Optimizer))

tpia_im_model = instantiate_mc_model(sf_math, FOTRUPowerModel, build_mc_tpia_L1; jump_model = Lower(bi_level_jump_model)) # ACPUPowerModel, OR LinDist3FlowPowerModel, OR FBSUBFPowerModel, OR FOTPUPowerModel, OR FOTRUPowerModel
tpia_jump_model = tpia_im_model.model
tpia_variables = all_variables(tpia_jump_model)
tpia_constraints = Vector{Any}()
for (function_type, set_type) in list_of_constraint_types(tpia_jump_model)
    append!(tpia_constraints, all_constraints(tpia_jump_model, function_type, set_type))
end

# run(`scp arnav.gautam@hammer.mines.edu:/home/arnav.gautam/model/ARNAV.nl ./test/data/`)
# jaime_model_NLP = AmplModel("test/data/ARNAV.nl")

# jaime_model_MOI = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_NL)
# MOI.read_from_file(jaime_model_MOI, "test/data/ARNAV.nl")

# INSTANTIATE THE JuMP MODEL FROM JAIME'S OPTIMIZATION FORMULATION

jaime_jump_model = JuMP.read_from_file("/home/arnav.gautam/model/ARNAV.nl")
jaime_variables = all_variables(jaime_jump_model)
jaime_constraints = Vector{Any}()
for (function_type, set_type) in list_of_constraint_types(jaime_jump_model)
    append!(jaime_constraints, all_constraints(jaime_jump_model, function_type, set_type))
end


# set_optimizer_attribute(bi_level_jump_model, "NonConvex", 2) # If I have Nonconvex/nonconcave/nonpsd objective/constraint error in a MIP solver like Gurobi.

@variable(Lower(bi_level_jump_model), tpia_variables[1:length(tpia_variables)])

for upper_var in jaime_variables
    base_name_value = name(upper_var)
    @variable(Upper(bi_level_jump_model), name(upper_var))
end

@variable(Upper(bi_level_jump_model), jaime_variables[1:length(jaime_variables)])

# set_optimizer(bi_level_jump_model, )

# for (function_type, set_type) in tpia_constrai

# Two ways to obtain convergence
# sf_tpia_solution = solve_mc_tpia_L1(sf_math, ACPUPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>200, "tol"=>100.0, "dual_inf_tol"=>100.0, "constr_viol_tol"=>1.0, "compl_inf_tol"=>0.01))
# sf_tpia_solution = solve_mc_tpia_L1(sf_math, ACPUPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>200, "acceptable_tol"=>100.0, "acceptable_dual_inf_tol"=>10000.0, "acceptable_constr_viol_tol"=>1.0, "acceptable_compl_inf_tol"=>0.1, "acceptable_obj_change_tol"=>10.0))

# Solve the bilevel optimization problem
# JuMP.optimize!(bi_level_jump_model)