"Three-Phase Infeasibility Analysis Problem"
function solve_mc_tpia(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_tpia; kwargs...)
end

"Three-Phase Infeasibility Analysis Problem with L1 Norm"
function solve_mc_tpia_L1(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_tpia_L1; kwargs...)
end

"Constructor for Three-Phase Infeasibility Analysis Problem"
function build_mc_tpia_L1(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage(pm; bounded=false)
    variable_mc_branch_power(pm; bounded=false)
    variable_mc_switch_power(pm; bounded=false)
    variable_mc_transformer_power(pm; bounded=false)
    variable_mc_generator_power(pm; bounded=false)
    @infiltrate
    variable_mc_load_power(pm; bounded=false)
    variable_mc_storage_power(pm; bounded=false)

    constraint_mc_model_voltage(pm)

    variable_mc_slack_bus_power_equity_weight(pm)
    variable_mc_slack_bus_power_L1(pm)

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_power_balance_slack_L1(pm, i)

        # PV Bus Constraints
        if (length(ref(pm, :bus_gens, i)) > 0 || length(ref(pm, :bus_storages, i)) > 0) && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
            for j in ref(pm, :bus_storages, i)
                constraint_mc_storage_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_nl(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_slack_bus_power_L1(pm)
end


"Constructor for Three-Phase Infeasibility Analysis Problem"
function build_mc_tpia(pm::AbstractUnbalancedPowerModel)
    # @infiltrate
    variable_mc_bus_voltage(pm; bounded=false)
    variable_mc_branch_power(pm; bounded=false)
    variable_mc_switch_power(pm; bounded=false)
    variable_mc_transformer_power(pm; bounded=false)
    variable_mc_generator_power(pm; bounded=false)
    variable_mc_load_power(pm; bounded=false)
    variable_mc_storage_power(pm; bounded=false)

    constraint_mc_model_voltage(pm)

    variable_mc_slack_bus_power_equity_weight(pm)
    variable_mc_slack_bus_power(pm)

    for (i,bus) in ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3

        constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_power_balance_slack(pm, i)
        # constraint_mc_slack_bus_power_equity_weight(pm, i)

        # PV Bus Constraints
        if (length(ref(pm, :bus_gens, i)) > 0 || length(ref(pm, :bus_storages, i)) > 0) && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
            for j in ref(pm, :bus_storages, i)
                constraint_mc_storage_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :storage)
        constraint_storage_state(pm, i)
        constraint_storage_complementarity_nl(pm, i)
        constraint_mc_storage_losses(pm, i)
        constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i)
    end

    objective_mc_min_slack_bus_power(pm)
end
