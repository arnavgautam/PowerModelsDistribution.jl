"Three-Phase Multi-Period Infeasibility Analysis Problem"
function solve_mc_tmpia(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_tmpia; kwargs...)
end

"Constructor for Three-Phase Multi-Period Infeasibility Analysis Problem"
function build_mc_tmpia(pm::AbstractUnbalancedPowerModel)
    
    # This variable should have one instance used in all constraints involving a cap on slack power
    base_nw = sort!(collect(keys(nws(pm))))[1]
    variable_mc_slack_bus_power_cap(pm; nw=base_nw)

    for (n, network) in nws(pm)
        println(string("network is ", n))
        variable_mc_bus_voltage(pm; nw=n, bounded=false)
        variable_mc_branch_power(pm; nw=n, bounded=false)
        variable_mc_switch_power(pm; nw=n, bounded=false)
        variable_mc_transformer_power(pm; nw=n, bounded=false)
        variable_mc_generator_power(pm; nw=n, bounded=false)
        variable_mc_load_power(pm; nw=n)
        variable_mc_storage_power(pm; nw=n, bounded=false)

        constraint_mc_model_voltage(pm; nw=n)

        variable_mc_slack_bus_power_L1(pm; nw=n)
        # variable_mc_slack_bus_power_equity_weight(pm) # Add back when incorporating as c_n variable
        
        for (i,bus) in ref(pm, :ref_buses; nw=n)
            @assert bus["bus_type"] == 3

            constraint_mc_theta_ref(pm, i; nw=n)
            constraint_mc_voltage_magnitude_only(pm, i; nw=n)
        end

        # gens should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :gen; nw=n)
            constraint_mc_generator_power(pm, id; nw=n)
        end

        # loads should be constrained before KCL, or Pd/Qd undefined
        for id in ids(pm, :load; nw=n)
            constraint_mc_load_power(pm, id; nw=n)
        end

        for (i,bus) in ref(pm, :bus; nw=n)
            constraint_mc_power_balance_slack_L1(pm, i; nw=n)
            constraint_mc_power_balance_slack_cap(pm, i; nw=n)

            # PV Bus Constraints
            if (length(ref(pm, :bus_gens, i; nw=n)) > 0 || length(ref(pm, :bus_storages, i; nw=n)) > 0) && !(i in ids(pm,:ref_buses; nw=n))
                # this assumes inactive generators are filtered out of bus_gens
                @assert bus["bus_type"] == 2

                constraint_mc_voltage_magnitude_only(pm, i; nw=n)
                for j in ref(pm, :bus_gens, i; nw=n)
                    constraint_mc_gen_power_setpoint_real(pm, j; nw=n)
                end
                for j in ref(pm, :bus_storages, i; nw=n)
                    constraint_mc_storage_power_setpoint_real(pm, j; nw=n)
                end
            end
        end

        for i in ids(pm, :storage; nw=n)
            constraint_storage_state(pm, i; nw=n)
            constraint_storage_complementarity_nl(pm, i; nw=n)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
        end

        for i in ids(pm, :branch; nw=n)
            constraint_mc_ohms_yt_from(pm, i; nw=n)
            constraint_mc_ohms_yt_to(pm, i; nw=n)
        end

        for i in ids(pm, :switch; nw=n)
            constraint_mc_switch_state(pm, i; nw=n)
        end

        for i in ids(pm, :transformer; nw=n)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
    end

    objective_mc_min_slack_bus_power_L1(pm)
end


####


"Three-Phase Multi-Period Infeasibility Analysis Problem"
function solve_mc_tmpia_single_nw(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_tmpia_single_nw; kwargs...)
end

"Constructor for Three-Phase Multi-Period Infeasibility Analysis Problem"
function build_mc_tmpia_single_nw(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage(pm; bounded=false)
    variable_mc_branch_power(pm; bounded=false)
    variable_mc_switch_power(pm; bounded=false)
    variable_mc_transformer_power(pm; bounded=false)
    variable_mc_generator_power(pm; bounded=false)
    variable_mc_load_power(pm)
    variable_mc_storage_power(pm; bounded=false)

    constraint_mc_model_voltage(pm)

    variable_mc_slack_bus_power_L1(pm)
    variable_mc_slack_bus_power_cap(pm)
    # variable_mc_slack_bus_power_equity_weight(pm) # Add back when incorporating as c_n variable
    
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
        constraint_mc_power_balance_slack_cap(pm, i)

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