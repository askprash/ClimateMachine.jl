using MPI
using ClimateMachine
using Logging
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.MPIStateArrays
using ClimateMachine.ODESolvers
using LinearAlgebra
using Printf
using Test
import ClimateMachine.VTK: writevtk, writepvtu
import ClimateMachine.GenericCallbacks:
    EveryXWallTimeSeconds, EveryXSimulationSteps

if !@isdefined integration_testing
    const integration_testing = parse(
        Bool,
        lowercase(get(ENV, "JULIA_CLIMA_INTEGRATION_TESTING", "false")),
    )
end

const output = parse(Bool, lowercase(get(ENV, "JULIA_CLIMA_OUTPUT", "false")))

include("advection_diffusion_model.jl")

struct Pseudo1D{n1, n2, α, β, μ, δ} <: AdvectionDiffusionProblem end

function init_velocity_diffusion!(
    ::Pseudo1D{n1, n2, α, β},
    aux::Vars,
    geom::LocalGeometry,
) where {n1, n2, α, β}
    # Direction of flow is n1 (resp n2) with magnitude α
    aux.advection.u = hcat(α * n1, α * n2)

    # diffusion of strength β in the n1 and n2 directions
    aux.diffusion.D = hcat(β * n1 * n1', β * n2 * n2')
end

function initial_condition!(
    ::Pseudo1D{n1, n2, α, β, μ, δ},
    state,
    aux,
    localgeo,
    t,
) where {n1, n2, α, β, μ, δ}
    ξn1 = dot(n1, localgeo.coord)
    ξn2 = dot(n2, localgeo.coord)
    ρ1 = exp(-(ξn1 - μ - α * t)^2 / (4 * β * (δ + t))) / sqrt(1 + t / δ)
    ρ2 = exp(-(ξn2 - μ - α * t)^2 / (4 * β * (δ + t))) / sqrt(1 + t / δ)
    state.ρ = (ρ1, ρ2)
end

Dirichlet_data!(P::Pseudo1D, x...) = initial_condition!(P, x...)

function Neumann_data!(
    ::Pseudo1D{n1, n2, α, β, μ, δ},
    ∇state,
    aux,
    x,
    t,
) where {n1, n2, α, β, μ, δ}
    ξn1 = dot(n1, x)
    ξn2 = dot(n2, x)
    ∇ρ1 =
        -(
            2n1 * (ξn1 - μ - α * t) / (4 * β * (δ + t)) *
            exp(-(ξn1 - μ - α * t)^2 / (4 * β * (δ + t))) / sqrt(1 + t / δ)
        )
    ∇ρ2 =
        -(
            2n2 * (ξn2 - μ - α * t) / (4 * β * (δ + t)) *
            exp(-(ξn2 - μ - α * t)^2 / (4 * β * (δ + t))) / sqrt(1 + t / δ)
        )
    ∇state.ρ = hcat(∇ρ1, ∇ρ2)
end

function do_output(mpicomm, vtkdir, vtkstep, dgfvm, Q, Qe, model, testname)
    ## name of the file that this MPI rank will write
    filename = @sprintf(
        "%s/%s_mpirank%04d_step%04d",
        vtkdir,
        testname,
        MPI.Comm_rank(mpicomm),
        vtkstep
    )

    statenames = flattenednames(vars_state(model, Prognostic(), eltype(Q)))
    exactnames = statenames .* "_exact"

    writevtk(filename, Q, dgfvm, statenames, Qe, exactnames)

    ## generate the pvtu file for these vtk files
    if MPI.Comm_rank(mpicomm) == 0
        ## name of the pvtu file
        pvtuprefix = @sprintf("%s/%s_step%04d", vtkdir, testname, vtkstep)

        ## name of each of the ranks vtk files
        prefixes = ntuple(MPI.Comm_size(mpicomm)) do i
            @sprintf("%s_mpirank%04d_step%04d", testname, i - 1, vtkstep)
        end

        writepvtu(
            pvtuprefix,
            prefixes,
            (statenames..., exactnames...),
            eltype(Q),
        )

        @info "Done writing VTK: $pvtuprefix"
    end
end


function test_run(mpicomm, dim, polynomialorders, level, ArrayType, FT, vtkdir)

    n_hd =
        dim == 2 ? SVector{3, FT}(1, 0, 0) :
        SVector{3, FT}(1 / sqrt(2), 1 / sqrt(2), 0)

    n_vd = dim == 2 ? SVector{3, FT}(0, 1, 0) : SVector{3, FT}(0, 0, 1)

    α = FT(1)
    β = FT(1 // 100)
    μ = FT(-1 // 2)
    δ = FT(1 // 10)

    # Grid/topology information
    base_num_elem = 4
    Ne = 2^(level - 1) * base_num_elem
    brickrange = ntuple(j -> range(FT(-1); length = Ne + 1, stop = 1), dim)
    periodicity = ntuple(j -> false, dim)
    bc = ntuple(j -> (1, 2), dim)

    topl = StackedBrickTopology(
        mpicomm,
        brickrange;
        periodicity = periodicity,
        boundary = bc,
    )

    dt = (α / 4) / (Ne * maximum(polynomialorders)^2)
    timeend = 1
    outputtime = timeend / 10
    @info "time step" dt

    @info @sprintf """Test parameters:
    ArrayType                   = %s
    FloatType                   = %s
    Dimension                   = %s
    Horizontal polynomial order = %s
    Vertical polynomial order   = %s
      """ ArrayType FT dim polynomialorders[1] polynomialorders[end]

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = ArrayType,
        polynomialorder = polynomialorders,
    )

    # Model being tested
    model = AdvectionDiffusion{dim}(
        Pseudo1D{n_hd, n_vd, α, β, μ, δ}(),
        num_equations = 2,
    )

    # Main DG discretization
    dgfvm = DGFVMModel(
        model,
        grid,
        RusanovNumericalFlux(),
        CentralNumericalFluxSecondOrder(),
        CentralNumericalFluxGradient();
        direction = EveryDirection(),
    )

    # Initialize all relevant state arrays and create solvers
    Q = init_ode_state(dgfvm, FT(0))

    eng0 = norm(Q, dims = (1, 3))
    @info @sprintf """Starting
    norm(Q₀) = %.16e""" eng0[1]

    solver = LSRK54CarpenterKennedy(dgfvm, Q; dt = dt, t0 = 0)

    # Set up the information callback
    starttime = Ref(Dates.now())
    cbinfo = EveryXWallTimeSeconds(60, mpicomm) do (s = false)
        if s
            starttime[] = Dates.now()
        else
            energy = norm(Q)
            @info @sprintf(
                """Update
                simtime = %.16e
                runtime = %s
                norm(Q) = %.16e""",
                gettime(solver),
                Dates.format(
                    convert(Dates.DateTime, Dates.now() - starttime[]),
                    Dates.dateformat"HH:MM:SS",
                ),
                energy
            )
        end
    end
    callbacks = (cbinfo,)
    if ~isnothing(vtkdir)
        # create vtk dir
        mkpath(vtkdir)

        vtkstep = 0
        # output initial step
        do_output(
            mpicomm,
            vtkdir,
            vtkstep,
            dgfvm,
            Q,
            Q,
            model,
            "advection_diffusion",
        )

        # setup the output callback
        cbvtk = EveryXSimulationSteps(floor(outputtime / dt)) do
            vtkstep += 1
            Qe = init_ode_state(dgfvm, gettime(solver))
            do_output(
                mpicomm,
                vtkdir,
                vtkstep,
                dgfvm,
                Q,
                Qe,
                model,
                "advection_diffusion",
            )
        end
        callbacks = (callbacks..., cbvtk)
    end

    solve!(Q, solver; timeend = timeend, callbacks = callbacks)

    # Reference solution
    engf = norm(Q, dims = (1, 3))
    Q_ref = init_ode_state(dgfvm, FT(timeend))

    engfe = norm(Q_ref, dims = (1, 3))
    errf = norm(Q_ref .- Q, dims = (1, 3))

    metrics = @. (engf, engf / eng0, engf - eng0, errf, errf / engfe)

    @info @sprintf """Finished
    Horizontal field:
      norm(Q)                 = %.16e
      norm(Q) / norm(Q₀)      = %.16e
      norm(Q) - norm(Q₀)      = %.16e
      norm(Q - Qe)            = %.16e
      norm(Q - Qe) / norm(Qe) = %.16e
    Vertical field:
      norm(Q)                 = %.16e
      norm(Q) / norm(Q₀)      = %.16e
      norm(Q) - norm(Q₀)      = %.16e
      norm(Q - Qe)            = %.16e
      norm(Q - Qe) / norm(Qe) = %.16e
      """ first.(metrics)... last.(metrics)...

    return errf
end

"""
    main()

Run this test problem
"""
function main()

    ClimateMachine.init()
    ArrayType = ClimateMachine.array_type()
    mpicomm = MPI.COMM_WORLD

    # Dictionary keys: dim, level, polynomial order, FT, and direction
    expected_result = Dict()

    @testset "Variable degree DG: advection diffusion model" begin
        for FT in (Float64,) #(Float32, Float64)
            numlevels =
                integration_testing ||
                ClimateMachine.Settings.integration_testing ?
                (FT == Float64 ? 4 : 3) : 1
            numlevels = 4
            for dim in 2:3
                polynomialorders = (4, 0)
                result = Dict()
                for level in 1:numlevels
                    vtkdir =
                        true || output ?
                        "vtk_advection" *
                        "_poly$(polynomialorders)" *
                        "_dim$(dim)_$(ArrayType)_$(FT)" *
                        "_level$(level)" :
                        nothing
                    result[level] = test_run(
                        mpicomm,
                        dim,
                        polynomialorders,
                        level,
                        ArrayType,
                        FT,
                        vtkdir,
                    )
                    #JK horiz_poly = polynomialorders[1]
                    #JK vert_poly = polynomialorders[2]
                    #JK @test result[level][1] ≈ FT(expected_result[
                    #JK     dim,
                    #JK     level,
                    #JK     horiz_poly,
                    #JK     FT,
                    #JK     HorizontalDirection,
                    #JK ])
                    #JK @test result[level][2] ≈ FT(expected_result[
                    #JK     dim,
                    #JK     level,
                    #JK     vert_poly,
                    #JK     FT,
                    #JK     VerticalDirection,
                    #JK ])
                end
                @info begin
                    msg = ""
                    for l in 1:(numlevels - 1)
                        rate = @. log2(result[l]) - log2(result[l + 1])
                        msg *= @sprintf(
                            "\n  rates for level %d Horizontal = %e",
                            l,
                            rate[1]
                        )
                        msg *= @sprintf(", Vertical = %e\n", rate[2])
                    end
                    msg
                end
            end
        end
    end
end

main()
