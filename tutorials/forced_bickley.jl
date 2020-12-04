# # Forced Bickley jet

using CUDA
using Printf
using Plots
using Revise
using ClimateMachine

ClimateMachine.init()

ClimateMachine.Settings.array_type = Array

using ClimateMachine.Ocean
using ClimateMachine.Ocean.Domains
using ClimateMachine.Ocean.Fields

using ClimateMachine.GenericCallbacks: EveryXSimulationTime
using ClimateMachine.GenericCallbacks: EveryXSimulationSteps
using ClimateMachine.Ocean: current_step, Δt, current_time
using ClimateMachine.Ocean: JLD2Writer, OutputTimeSeries, write!
using CLIMAParameters: AbstractEarthParameterSet, Planet

# Domain

domain = RectangularDomain(
    Ne = (64, 64, 1), Np = 4,
    x = (-2π, 2π), y = (-2π, 2π), z = (0, 1),
    periodicity = (true, false, false),
)

# Non-dimensional parameters:

struct NonDimensionalParameters <: AbstractEarthParameterSet end
Planet.grav(::NonDimensionalParameters) = 10
g = Planet.grav(NonDimensionalParameters())

f = 1e-3 # Coriolis parameter
ϵ = 0.2 # Perturbation amplitude
ℓ = 0.5 # Perturbation width
k = 2.0 # Perturbation wavenumber
const a = f / g # Surface displacement amplitude

# Jet
Ψ(y) = - tanh(y)
U(y) = sech(y)^2

# Tracer
const L = domain.L.x
T(y) = sin(2π * y / L)

# Perturbations
ψ̃(x, y) = exp(-y^2 / 2ℓ^2) * cos(k * x) * cos(k * y)

# ũ = - ∂_y ψ̃ = (k * tan(k * y) + y / ℓ^2) * ψ̃
ũ(x, y) = (k * tan(k * y) + y / ℓ^2) * ψ̃(x, y)

# ṽ = + ∂_x ψ̃ = - k * tan(k * x) * ψ̃
ṽ(x, y) = - k * tan(k * x) * ψ̃(x, y)

# Initial conditions: Jet/tracer + perturbations
uᵢ(x, y, z) = U(y) + ϵ * ũ(x, y)
vᵢ(x, y, z) = ϵ * ṽ(x, y)
θᵢ(x, y, z) = T(y)
ηᵢ(x, y, z) = a * (Ψ(y) + ϵ * ψ̃(x, y))

initial_conditions = InitialConditions(u=uᵢ, v=vᵢ, θ=θᵢ, η=ηᵢ)

# Forcing: relax to large-scale state
const τ = 100.0
relax_u(y, t, u, v, w, η, θ) = 1 / τ * (U(y) - u)
relax_η(y, t, u, v, w, η, θ) = 1 / τ * (a * Ψ(y) - η)
relax_θ(y, t, u, v, w, η, θ) = 1 / τ * (T(y) - θ)

forcing = Ocean.Forcing(u=relax_u, η=relax_η, θ=relax_θ)

model = Ocean.HydrostaticBoussinesqSuperModel(
    domain = domain,
    time_step = 0.002,
    initial_conditions = initial_conditions,
    array_type = CuArray,
    parameters = NonDimensionalParameters(),
    turbulence_closure = (νʰ = 1e-4, κʰ = 1e-4,
                          νᶻ = 1e-2, κᶻ = 1e-2),
    forcing = forcing,
    rusanov_wave_speeds = (cʰ = sqrt(g * domain.L.z), cᶻ = 1e-2),
    coriolis = (f₀ = f, β = 0),
    buoyancy = (αᵀ = 0,),
    boundary_tags = ((0, 0), (1, 1), (1, 2)),
    boundary_conditions = (
        OceanBC(Impenetrable(FreeSlip()), Insulating()),
        OceanBC(Penetrable(FreeSlip()), Insulating()),
    ),
)

# We prepare a callback that periodically fetches the horizontal velocity and
# tracer concentration for later animation,

realdata = convert(Array, model.state.realdata)
u = SpectralElementField(domain, model.grid, view(realdata, :, 1, :))

volume = assemble(u)
x = volume.x[:, 1, 1]
y = volume.y[1, :, 1]

writer = JLD2Writer(model, filepath="test.jld2", overwrite_existing=true)

start_time = time_ns()

data_fetcher = EveryXSimulationTime(1) do
    write!(writer)

    realdata = convert(Array, model.state.realdata)
    u = SpectralElementField(domain, model.grid, view(realdata, :, 1, :))

    # Print a helpful message
    step = @sprintf("Step: %d", current_step(model))
    time = @sprintf("time: %.2f", current_time(model))
    max_u = @sprintf("max|u|: %.6f", maximum(abs, u))

    elapsed = (time_ns() - start_time) * 1e-9
    wall_time = @sprintf("elapsed wall time: %.2f min", elapsed / 60)  

    isnan(maximum(abs, u)) && error("NaNs.") 

    @info "$step, $time, $max_u, $wall_time"
end

# and then run the simulation.

model.solver_configuration.timeend = 200

try
    result = ClimateMachine.invoke!(
        model.solver_configuration;
        user_callbacks = [data_fetcher],
    )
catch err
    @warn "Simulation ended prematurely because $(sprint(showerror, err))"
end

# Finally, we make an animation of the evolving shear instability.

timeserieses = OutputTimeSeries.(keys(model.fields), writer.filepath)

u_timeseries = OutputTimeSeries(:u, writer.filepath)
v_timeseries = OutputTimeSeries(:v, writer.filepath)
η_timeseries = OutputTimeSeries(:η, writer.filepath)
θ_timeseries = OutputTimeSeries(:θ, writer.filepath)

animation = @animate for i = 1:length(u_timeseries)

    local u

    @info "Plotting frame $i of $(length(u_timeseries))..."

    kwargs = (xlim = domain.x, ylim = domain.y, linewidth = 0, aspectratio = 1)

    u = assemble(u_timeseries[i]).data[:, :, 1]
    v = assemble(v_timeseries[i]).data[:, :, 1]
    η = assemble(η_timeseries[i]).data[:, :, 1]
    θ = assemble(θ_timeseries[i]).data[:, :, 1]

    ulim = 1
    θlim = 1

    ηmin = minimum(η)
    ηmax = maximum(η)

    ulevels = range(-ulim, ulim, length=31)
    ηlevels = range(ηmin, ηmax, length=31)
    θlevels = range(-θlim, θlim, length=31)

    u_plot = heatmap(x, y, clamp.(u, -ulim, ulim)';
                     levels = ulevels, clim = (-ulim, ulim), color = :balance, kwargs...)

    v_plot = heatmap(x, y, clamp.(v, -ulim, ulim)';
                     levels = ulevels, clim = (-ulim, ulim), color = :balance, kwargs...)

    η_plot = heatmap(x, y, clamp.(η, ηmin, ηmax)';
                     levels = ηlevels, clim = (ηmin, ηmax), color = :balance, kwargs...)

    θ_plot = heatmap(x, y, clamp.(θ, -θlim, θlim)';
                     levels = θlevels, clim = (-θlim, θlim), color = :thermal, kwargs...)

    u_title = @sprintf("u at t = %.2f", u_timeseries.times[i])
    v_title = @sprintf("v at t = %.2f", u_timeseries.times[i])
    θ_title = @sprintf("θ at t = %.2f", u_timeseries.times[i])
    η_title = @sprintf("η at t = %.2f", u_timeseries.times[i])

    plot(u_plot, v_plot, η_plot, θ_plot, title = [u_title v_title η_title θ_title],
         layout = (2, 2), size = (1200, 1000))
end

gif(animation, "forced_bickley.gif", fps = 8)
