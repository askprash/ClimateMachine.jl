#!/usr/bin/env julia --project
using ClimateMachine

using ClimateMachine.Atmos
using ClimateMachine.Orientations
using ClimateMachine.ConfigTypes
using ClimateMachine.Diagnostics
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.GenericCallbacks
using ClimateMachine.ODESolvers
using ClimateMachine.Mesh.Filters
using ClimateMachine.TemperatureProfiles
using ClimateMachine.Thermodynamics
using ClimateMachine.TurbulenceClosures
using ClimateMachine.VariableTemplates
using ClimateMachine.BalanceLaws:
    AbstractStateType, Auxiliary, UpwardIntegrals, DownwardIntegrals

using ArgParse
using Distributions
using Random
using StaticArrays
using Test
using DocStringExtensions
using LinearAlgebra
using Printf

using CLIMAParameters
using CLIMAParameters.Planet: planet_radius, R_d, Omega, cp_d, MSLP, grav, LH_v0, cv_d

using CLIMAParameters.Atmos.Microphysics

struct LiquidParameterSet <: AbstractLiquidParameterSet end
struct IceParameterSet <: AbstractIceParameterSet end
struct RainParameterSet <: AbstractRainParameterSet end
struct SnowParameterSet <: AbstractSnowParameterSet end

struct MicropysicsParameterSet{L, I, R, S} <: AbstractMicrophysicsParameterSet
    liq::L
    ice::I
    rai::R
    sno::S
end

struct EarthParameterSet{M} <: AbstractEarthParameterSet
    microphys::M
end

microphys = MicropysicsParameterSet(
    LiquidParameterSet(),
    IceParameterSet(),
    RainParameterSet(),
    SnowParameterSet(),
)

const param_set = EarthParameterSet(microphys)


import ClimateMachine.BalanceLaws:
    vars_state,
    indefinite_stack_integral!,
    reverse_indefinite_stack_integral!,
    integral_load_auxiliary_state!,
    integral_set_auxiliary_state!,
    reverse_integral_load_auxiliary_state!,
    reverse_integral_set_auxiliary_state!

import ClimateMachine.BalanceLaws: boundary_state!
import ClimateMachine.Atmos: flux_second_order!

"""
  Initial Condition for DYCOMS_RF01 LES

## References
 - [Stevens2005](@cite)
"""
function init_dycoms!(problem, bl, state, aux, localgeo, t)
    FT = eltype(state)

    (x, y, z) = localgeo.coord

    z = altitude(bl, aux)

    # These constants are those used by Stevens et al. (2005)
    q0 = FT(0.021)
    qtrop = 1e-11
    zq1 = FT(3000)
    zq2 = FT(8000)
    theta_0 = FT(1)
    sigma_r = FT(20000)
    sigma_z = FT(2000)
    r_b = FT(40000)
    z_b = FT(5000)
    q_pt_sfc = PhasePartition(qref)
    Rm_sfc = gas_constant_air(bl.param_set, q_pt_sfc)
    Ts0 = FT(302.15)
    constTv = FT(0.608)
    _MSLP = FT(MSLP(bl.param_set))
    _grav = FT(grav(bl.param_set))
    rp = FT(282000)
    cen_lat = FT(0)
    cen_lon = FT(0)
    convert = FT(180) / pi
    _Omega::FT = Omega(m.param_set)
    _R_d::FT = R_d(m.param_set)
    _cp::FT = cp_d(m.param_set)
    _cv::FT = cv_d(m.param_set)
    γ::FT = _cp / _cv
    _PR = planet_radius(m.param_set)
    ztrop = FT(15000)
    exppr = FT(1.5)
    zp = FT(7000)
    exppz = FT(2)
    exponent = _R_d * γ / _grav
    T0 = Ts0 * (FT(1) + constTv * q0)
    Ttrop = Ts0 - γ * ztrop
    ptrop = _MSLP * (Ttrop / T0)^(FT(1) / exponent)
    lat = y
    lon = x
    
    f = FT(2) * _Omega * sin(cen_lat / convert)
    gr = _PR * acos(sin(cen_lat / convert) * sin(lat) + (cos(cen_lat / convert) * cos(lat) * cos(lon - cen_lon / convert))) 
    # if cartesian
    gr = sqrt((cen_lat - lat)^2 + (cen_lon - lon)^2)
    f = FT(5e-5)
    ps = _MSLP - dp * exp(-(gr / rp)^exppr)
    height = z
    if (height > ztrop)
      p = ptrop * exp(-(_grav * (height - ztrop)) / (_R_d * Ttrop))
      pb = p
    else 
      p = (_MSLP - dp * exp(-(gr / rp)^exppr) * exp(-(height / zp)^exppz)) * ((T0 - γ * height) / T0)^(1/exponent)
      pb = _MSLP * ((T0 - γ * height) / T0)^(1/exponent)
    end
    d1 = sin(cen_lat / convert) * cos(lat) - cos(cen_lat / convert) * sin(lat) * cos(lon - cen_lon / convert)
    d2 = cos(cen_lat / convert) * sin(lon - cen_lon / convert)
    d = max(1e-25, sqrt(d1^2 + d2^2))
    ufac = d1 / d
    vfac = d2 / d
    if not sphere
      angle = atan(lon - cen_lon, lat - cen_lat)
      theta = pi / 2 + angle
      radial_decay = exp(-height * height/(2 * 5823 * 5823)) * exp(-(gr/200000)^6)
      ufac = cos(theta) * radial_decay
      vfac = sin(theta) * radial_decay
    end
    if (height > ztrop)
      us = FT(0)
      vs = FT(0)
    else 
      vs = vfac * (-f * gr / 2 + sqrt((f * gr / 2)^2 - exppr * (gr / rp)^exppr * _R_d * (T0 - γ * height) / (exppz * height * R_d * (T0 - γ * height) / (_grav * zp^exppz) + (FT(1) - _MSLP / dp * exp((gr/rp)^exppr) * exp((height / zp)^exppz)))))
      us = ufac * (-f * gr / 2 + sqrt((f * gr / 2)^2 - exppr * (gr / rp)^exppr * _R_d * (T0 - γ * height) / (exppz * height * R_d * (T0 - γ * height) / (_grav * zp^exppz) + (FT(1) - _MSLP / dp * exp((gr/rp)^exppr) * exp((height / zp)^exppz)))))
    end
    u = FT(us)
    v = FT(vs)
    w = FT(0) 
    if sphere
      u = -us * sin(lon) - vs * sin(lat) * cos(lon)
      v = us * cos(lon) - vs * sin(lat) * cos(lon)
      w = vs * cos(lat)
    end
    if (height > ztrop)
       q = qtrop
    else 
       q = q0*exp(-height/zq1) * exp(-(height/zq2)^exppz)
    end
    qb = q
    if (height > ztrop)
      T = Ttrop
      Tb = T
    else
      T = (T0 - γ * height) / (FT(1) + constTv * q) / (FT(1) + exppz * _R_d * (T0 - γ * height) * height / (_grav * zp * exppz * (FT(1) - _MSLP / dp * exp((gr/rp)^exppr) * exp((height/zp)^exppz))))
      Tb = (T0 - γ * height)/(FT(1) + constTv * q)
    end
    wavenumber = 3
    angle = atan(lat - cen_lat, lon - cen_lon)
    wave = FT(exp(complex(0,1) * (wavenumber * (pi / 2 - angle))))
    pert = wave * theta_0 * exp(-((gr - r_b) / sigma_r)^2 - ((height - z_b) / sigma_z)^2)
    T = T + pert

    ts = PhaseEquil_pθq(bl.param_set, p, θ_liq, q_tot)
    ρ = air_density(ts)

    e_kin = FT(1 / 2) * FT((u^2 + v^2 + w^2))
    e_pot = gravitational_potential(bl.orientation, aux)
    E = ρ * total_energy(e_kin, e_pot, ts)

    state.ρ = ρ
    state.ρu = SVector(ρ * u, ρ * v, ρ * w)
    state.ρe = E

    state.moisture.ρq_tot = ρ * q_tot

    if bl.moisture isa NonEquilMoist
        q_init = PhasePartition(ts)
        state.moisture.ρq_liq = q_init.liq
        state.moisture.ρq_ice = q_init.ice
    end
    if bl.precipitation isa Rain
        state.precipitation.ρq_rai = FT(0)
    end

    return nothing
end

function config_dycoms(
    FT,
    N,
    resolution,
    xmax,
    ymax,
    zmax,
    moisture_model = "equilibrium",
    precipitation_model = "noprecipitation",
)
    # Reference state
    T_profile = DecayingTemperatureProfile{FT}(param_set)
    ref_state = HydrostaticState(T_profile)

    # Radiation model
    κ = FT(85)
    α_z = FT(1)
    z_i = FT(840)
    ρ_i = FT(1.13)

    D_subsidence = FT(3.75e-6)

    F_0 = FT(70)
    F_1 = FT(22)
    if moisture_model == "equilibrium"
        equilibrium_moisture_model = true
    else
        equilibrium_moisture_model = false
    end
    radiation = DYCOMSRadiation{FT}(
        κ,
        α_z,
        z_i,
        ρ_i,
        D_subsidence,
        F_0,
        F_1,
        equilibrium_moisture_model,
    )

    # Sources
    f_coriolis = FT(0.762e-4)
    u_geostrophic = FT(7.0)
    v_geostrophic = FT(-5.5)
    w_ref = FT(0)
    u_relaxation = SVector(u_geostrophic, v_geostrophic, w_ref)
    # Sponge
    c_sponge = 1
    # Rayleigh damping
    zsponge = FT(1000.0)
    rayleigh_sponge =
        RayleighSponge(FT, zmax, zsponge, c_sponge, u_relaxation, 2)
    # Geostrophic forcing
    geostrophic_forcing =
        GeostrophicForcing(FT, f_coriolis, u_geostrophic, v_geostrophic)

    # Boundary conditions
    # SGS Filter constants
    C_smag = FT(0.21) # 0.21 for stable testing, 0.18 in practice
    C_drag = FT(0.0011)
    LHF = FT(115)
    SHF = FT(15)
    moisture_flux = LHF / FT(LH_v0(param_set))

    source = (
        Gravity(),
        rayleigh_sponge,
        Subsidence(D_subsidence)...,
        geostrophic_forcing,
    )

    # moisture model and its sources
    if moisture_model == "equilibrium"
        moisture = EquilMoist{FT}(; maxiter = 4, tolerance = FT(1))
    elseif moisture_model == "nonequilibrium"
        source = (source..., CreateClouds()...)
        moisture = NonEquilMoist()
    else
        @warn @sprintf(
            """
%s: unrecognized moisture_model in source terms, using the defaults""",
            moisture_model,
        )
        moisture = EquilMoist{FT}(; maxiter = 4, tolerance = FT(1))
    end

    # precipitation model and its sources
    if precipitation_model == "noprecipitation"
        precipitation = NoPrecipitation()
    elseif precipitation_model == "rain"
        source = (source..., Rain_1M())
        precipitation = Rain()
    else
        @warn @sprintf(
            """
%s: unrecognized precipitation_model in source terms, using the defaults""",
            precipitation_model,
        )
        precipitation = NoPrecipitation()
    end

    problem = AtmosProblem(
        boundaryconditions = (
            AtmosBC(
                momentum = Impenetrable(DragLaw(
                    (state, aux, t, normPu) -> C_drag,
                )),
                energy = PrescribedEnergyFlux((state, aux, t) -> LHF + SHF),
                moisture = PrescribedMoistureFlux(
                    (state, aux, t) -> moisture_flux,
                ),
            ),
            AtmosBC(),
        ),
        init_state_prognostic = init_dycoms!,
    )

    model = AtmosModel{FT}(
        AtmosLESConfigType,
        param_set;
        problem = problem,
        ref_state = ref_state,
        turbulence = Vreman{FT}(C_smag),
        moisture = moisture,
        precipitation = precipitation,
        radiation = radiation,
        source = source,
    )

    ode_solver = ClimateMachine.ExplicitSolverType(
        solver_method = LSRK144NiegemannDiehlBusch,
    )

    config = ClimateMachine.AtmosLESConfiguration(
        "DYCOMS",
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        param_set,
        init_dycoms!,
        solver_type = ode_solver,
        model = model,
    )
    return config
end

function config_diagnostics(driver_config)
    interval = "10000steps"
    dgngrp = setup_atmos_default_diagnostics(
        AtmosLESConfigType(),
        interval,
        driver_config.name,
    )
    return ClimateMachine.DiagnosticsConfiguration([dgngrp])
end

function main()
    # add a command line argument to specify the kind of
    # moisture and precipitation model you want
    # TODO: this will move to the future namelist functionality
    dycoms_args = ArgParseSettings(autofix_names = true)
    add_arg_group!(dycoms_args, "DYCOMS")
    @add_arg_table! dycoms_args begin
        "--moisture-model"
        help = "specify cloud condensate model"
        metavar = "equilibrium|nonequilibrium"
        arg_type = String
        default = "equilibrium"
        "--precipitation-model"
        help = "specify precipitation model"
        metavar = "noprecipitation|rain"
        arg_type = String
        default = "noprecipitation"
    end

    cl_args =
        ClimateMachine.init(parse_clargs = true, custom_clargs = dycoms_args)
    moisture_model = cl_args["moisture_model"]
    precipitation_model = cl_args["precipitation_model"]

    FT = Float64

    # DG polynomial order
    N = 4

    # Domain resolution and size
    Δh = FT(40)
    Δv = FT(20)
    resolution = (Δh, Δh, Δv)

    xmax = FT(1000)
    ymax = FT(1000)
    zmax = FT(1500)

    t0 = FT(0)
    timeend = FT(100) #FT(4 * 60 * 60)
    Cmax = FT(1.7)     # use this for single-rate explicit LSRK144

    driver_config = config_dycoms(
        FT,
        N,
        resolution,
        xmax,
        ymax,
        zmax,
        moisture_model,
        precipitation_model,
    )
    solver_config = ClimateMachine.SolverConfiguration(
        t0,
        timeend,
        driver_config,
        init_on_cpu = true,
        Courant_number = Cmax,
    )
    dgn_config = config_diagnostics(driver_config)

    if moisture_model == "equilibrium"
        filter_vars = ("moisture.ρq_tot",)
    elseif moisture_model == "nonequilibrium"
        filter_vars = ("moisture.ρq_tot", "moisture.ρq_liq", "moisture.ρq_ice")
    end
    if precipitation_model == "rain"
        filter_vars = (filter_vars..., "precipitation.ρq_rai")
    end

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do
        Filters.apply!(
            solver_config.Q,
            filter_vars,
            solver_config.dg.grid,
            TMARFilter(),
        )
        nothing
    end

    result = ClimateMachine.invoke!(
        solver_config;
        diagnostics_config = dgn_config,
        user_callbacks = (cbtmarfilter,),
        check_euclidean_distance = true,
    )
end

main()
