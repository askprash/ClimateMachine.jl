@horizontal_average(AtmosLESConfigType, u, "m s^-1", "x-velocity", "")
@horizontal_average(AtmosLESConfigType, v, "m s^-1", "y-velocity", "")
@horizontal_average(AtmosLESConfigType, w, "m s^-1", "z-velocity", "")

@horizontal_average(AtmosLESConfigType, uu)
@horizontal_average(AtmosLESConfigType, vv)
@horizontal_average(AtmosLESConfigType, ww)
@horizontal_average(AtmosLESConfigType, www)
@horizontal_average(AtmosLESConfigType, wu)
@horizontal_average(AtmosLESConfigType, wv)

@horizontal_average_impl(
    AtmosLESConfigType,
    u,
    v,
    w,
    uu,
    vv,
    ww,
    www,
    wu,
    wv,
) do (atmos::AtmosModel, states::States, curr_time)
    u = states.prognostic.ρu[1]
    v = states.prognostic.ρu[2]
    w = states.prognostic.ρu[3]
    uu = u^2 / states.prognostic.ρ
    vv = v^2 / states.prognostic.ρ
    ww = w^2 / states.prognostic.ρ
    www = w^3 / states.prognostic.ρ^2
    wu = w * u / states.prognostic.ρ
    wv = w * v / states.prognostic.ρ
    return u, v, w, uu, vv, ww, www, wu, wv
end

@horizontal_average(
    AtmosLESConfigType,
    avg_rho,
    "kg m^-3",
    "air density",
    "air_density"
)
@horizontal_average(
    AtmosLESConfigType,
    rho,
    "kg m^-3",
    "air density",
    "air_density"
)

@horizontal_average(AtmosLESConfigType, "wrho")

@horizontal_average_impl(
    AtmosLESConfigType,
    avg_rho,
    rho,
    wrho,
) do (atmos::AtmosModel, states::States, curr_time)
    avg_rho = states.prognostic.ρ
    rho = states.prognostic.ρ * states.prognostic.ρ
    wrho = states.prognostic.ρu[3] * states.prognostic.ρ
    return avg_rho, rho, wrho
end

@horizontal_average(
    AtmosLESConfigType,
    temp,
    "K",
    "air temperature",
    "air_temperature"
)
@horizontal_average(
    AtmosLESConfigType,
    pres,
    "Pa",
    "air pressure",
    "air_pressure"
)
@horizontal_average(
    AtmosLESConfigType,
    thd,
    "K",
    "dry potential temperature",
    "air_potential_temperature"
)
@horizontal_average(
    AtmosLESConfigType,
    et,
    "J kg^-1",
    "total specific energy",
    "specific_dry_energy_of_air"
)
@horizontal_average(
    AtmosLESConfigType,
    ei,
    "J kg^-1",
    "specific internal energy",
    "internal_energy"
)
@horizontal_average(
    AtmosLESConfigType,
    ht,
    "J kg^-1",
    "specific enthalpy based on total energy",
    ""
)
@horizontal_average(
    AtmosLESConfigType,
    hi,
    "J kg^-1",
    "specific enthalpy based on internal energy",
    "atmosphere_enthalpy_content"
)
@horizontal_average(
    AtmosLESConfigType,
    w_ht_sgs,
    "kg kg^-1 m s^-1",
    "vertical sgs flux of total specific enthalpy",
    ""
)

@horizontal_average(AtmosLESConfigType, eiei)
@horizontal_average(AtmosLESConfigType, wthd)
@horizontal_average(AtmosLESConfigType, wei)

@horizontal_average_impl(
    AtmosLESConfigType,
    temp,
    pres,
    thd,
    et,
    ei,
    ht,
    hi,
    w_ht_sgs,
    eiei,
    wthd,
    wei,
) do (atmos::AtmosModel, states::States, curr_time)
    e_tot = states.prognostic.ρe / states.prognostic.ρ
    ts = recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    air_temp = air_temperature(ts)
    air_pres = air_pressure(ts)
    θ_dry = dry_pottemp(ts)
    e_int = internal_energy(ts)
    h_tot = total_specific_enthalpy(ts, e_tot)
    h_int = specific_enthalpy(ts)

    temp = air_temp * states.prognostic.ρ
    pres = air_pres * states.prognostic.ρ
    thd = θ_dry * states.prognostic.ρ
    et = states.prognostic.ρe
    ei = e_int * states.prognostic.ρ
    ht = h_tot * states.prognostic.ρ
    hi = h_int * states.prognostic.ρ
    _, D_t, _ = turbulence_tensors(
        atmos,
        states.prognostic,
        states.gradient_flux,
        states.auxiliary,
        curr_time,
    )
    d_h_tot = -D_t .* states.gradient_flux.∇h_tot
    w_ht_sgs = d_h_tot[end] * states.prognostic.ρ
    eiei = ei * e_int
    wthd = states.prognostic.ρu[3] * θ_dry
    wei = states.prognostic.ρu[3] * e_int
    return temp, pres, thd, et, ei, ht, hi, w_ht_sgs, eiei, wthd, wei
end

@horizontal_average(
    AtmosLESConfigType,
    qt,
    "kg kg^-1",
    "mass fraction of total water in air (qv+ql+qi)",
    "mass_fraction_of_water_in_air"
)
@horizontal_average(
    AtmosLESConfigType,
    ql,
    "kg kg^-1",
    "mass fraction of liquid water in air",
    "mass_fraction_of_cloud_liquid_water_in_air"
)
@horizontal_average(
    AtmosLESConfigType,
    qi,
    "kg kg^-1",
    "mass fraction of ice in air",
    "mass_fraction_of_cloud_ice_in_air"
)
@horizontal_average(
    AtmosLESConfigType,
    qv,
    "kg kg^-1",
    "mass fraction of water vapor in air",
    "specific_humidity"
)
@horizontal_average(
    AtmosLESConfigType,
    thv,
    "K",
    "virtual potential temperature",
    "virtual_potential_temperature"
)
@horizontal_average(
    AtmosLESConfigType,
    thl,
    "K",
    "liquid-ice potential temperature",
    ""
)
@horizontal_average(
    AtmosLESConfigType,
    w_qt_sgs,
    "kg kg^-1 m s^-1",
    "vertical sgs flux of total specific humidity",
    ""
)

@horizontal_average(AtmosLESConfigType, qtqt)
@horizontal_average(AtmosLESConfigType, thlthl)
@horizontal_average(AtmosLESConfigType, wqt)
@horizontal_average(AtmosLESConfigType, wql)
@horizontal_average(AtmosLESConfigType, wqi)
@horizontal_average(AtmosLESConfigType, wqv)
@horizontal_average(AtmosLESConfigType, wthv)
@horizontal_average(AtmosLESConfigType, wthl)
@horizontal_average(AtmosLESConfigType, qtthl)
@horizontal_average(AtmosLESConfigType, qtei)

@horizontal_average_impl(
    AtmosLESConfigType,
    qt,
    ql,
    qi,
    qv,
    thv,
    thl,
    w_qt_sgs,
    qtqt,
    thlthl,
    wqt,
    wql,
    wqi,
    wqv,
    wthv,
    wthl,
    qtthl,
    qtei,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
)
    ts = recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    e_int = internal_energy(ts)
    q_liq = liquid_specific_humidity(ts)
    q_ice = ice_specific_humidity(ts)
    q_vap = vapor_specific_humidity(ts)
    θ_vir = virtual_pottemp(ts)
    θ_liq_ice = liquid_ice_pottemp(ts)
    has_condensate = has_condensate(ts)

    qt = states.prognostic.moisture.ρq_tot
    ql = q_liq * states.prognostic.ρ
    qi = q_ice * states.prognostic.ρ
    qv = q_vap * states.prognostic.ρ
    thv = θ_vir * states.prognostic.ρ
    thl = θ_liq_ice * states.prognostic.ρ
    _, D_t, _ = turbulence_tensors(
        atmos,
        states.prognostic,
        states.gradient_flux,
        states.auxiliary,
        curr_time,
    )
    d_q_tot = -D_t .* states.gradient_flux.moisture.∇q_tot
    w_qt_sgs = d_q_tot[end] * states.prognostic.ρ
    qtqt = qt * (states.prognostic.moisture.ρq_tot / states.prognostic.ρ)
    thlthl = thl * θ_liq_ice
    w = states.prognostic.ρu[3] / states.prognostic.ρ
    wqt = states.prognostic.ρu[3] * qt / states.prognostic.ρ
    wql = states.prognostic.ρu[3] * q_liq
    wqi = states.prognostic.ρu[3] * q_ice
    wqv = states.prognostic.ρu[3] * q_vap
    wthv = states.prognostic.ρu[3] * θ_vir
    wthl = states.prognostic.ρu[3] * θ_liq_ice
    qtthl = qt * θ_liq_ice
    qtei = qt * e_int
    return qt,
    ql,
    qi,
    qv,
    thv,
    thl,
    w_qt_sgs,
    qtqt,
    thlthl,
    wqt,
    wql,
    wqi,
    wqv,
    wthv,
    wthl,
    qtthl,
    qtei
end

#= TODO
    Variables["cld_frac"] = DiagnosticVariable(
        "cld_frac",
        diagnostic_var_attrib(
            "",
            "cloud fraction",
            "cloud_area_fraction_in_atmosphere_layer",
        ),
    )
    Variables["cld_cover"] = DiagnosticVariable(
        "cld_cover",
        diagnostic_var_attrib("", "cloud cover", "cloud_area_fraction"),
    )
    Variables["cld_top"] = DiagnosticVariable(
        "cld_top",
        diagnostic_var_attrib("m", "cloud top", "cloud_top_altitude"),
    )
    Variables["cld_base"] = DiagnosticVariable(
        "cld_base",
        diagnostic_var_attrib("m", "cloud base", "cloud_base_altitude"),
    )
    Variables["lwp"] = DiagnosticVariable(
        "lwp",
        diagnostic_var_attrib(
            "kg m^-2",
            "liquid water path",
            "atmosphere_mass_content_of_cloud_condensed_water",
        ),
    )
=#
