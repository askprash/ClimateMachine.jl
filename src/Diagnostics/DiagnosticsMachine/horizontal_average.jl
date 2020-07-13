"""
    HorizontalAverage

A horizontal reduction into a single vertical dimension.
"""
abstract type HorizontalAverage <: DiagnosticVar end
dv_HorizontalAverage(
    ::ClimateMachineConfigType,
    ::Union{HorizontalAverage},
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
) = nothing

dv_dimnames(::ClimateMachineConfigType, ::HorizontalAverage, ::Any) = ("z",)

macro horizontal_average(config_type, name)
    iex = generate_dv_interface(:HorizontalAverage, config_type, name)
    esc(MacroTools.prewalk(unblock, iex))
end

macro horizontal_average(config_type, name, units, long_name, standard_name)
    iex = generate_dv_interface(
        :HorizontalAverage,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )
    esc(MacroTools.prewalk(unblock, iex))
end

macro horizontal_average(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = quote
        $(generate_dv_interface(
            :HorizontalAverage,
            config_type,
            name,
            units,
            long_name,
            standard_name,
        ))
        $(generate_dv_function(:HorizontalAverage, config_type, [name], impl))
    end
    esc(MacroTools.prewalk(unblock, iex))
end

macro horizontal_averages(impl, config_type, names...)
    exprs = [
        generate_dv_interface(:HorizontalAverage, config_type, name)
        for name in names
    ]
    fex = generate_dv_function(:HorizontalAverage, config_type, names, impl)
    push!(exprs, fex)
    iex = quote
        $(exprs...)
    end
    esc(MacroTools.prewalk(unblock, iex))
end

macro horizontal_average_impl(impl, config_type, names...)
    iex = generate_dv_function(:HorizontalAverage, config_type, names, impl)
    esc(MacroTools.prewalk(unblock, iex))
end
