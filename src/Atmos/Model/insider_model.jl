# Insider Model allows computation of variables using model kernels using mini balance laws
using DocStringExtensions
using ..TemperatureProfiles
export ReferenceState, NoReferenceState, HydrostaticState
const TD = Thermodynamics
using CLIMAParameters.Planet: R_d, MSLP, cp_d, grav, T_surf_ref, T_min_ref

using ClimateMachine.BalanceLaws:
    AbstractStateType, Auxiliary, Gradient

import ClimateMachine.BalanceLaws:
    vars_state

"""
    VorticityModel

A mini balance law that is used to take the gradient of u and v to obtain vert

reference
pressure. The gradient is computed as ∇ ⋅(pI) and the calculation
uses the balance law interface to be numerically consistent with
the way this gradient is computed in the dynamics.
"""

struct VorticityModel{M} <: BalanceLaw
    atmos::M
end
  
function vars_state(m::VorticityModel, FT)
    @vars begin
        ω::SVector{3,FT}
    end
end
  
  vars_gradient(m::VorticityModel, FT) = @vars()
  vars_diffusive(m::VorticityModel, FT) = @vars()
  vars_aux(m::VorticityModel, FT) = vars_state(m.atmos, FT)


vars_state(::VorticityModel, ::Auxiliary, T) = @vars(ω::T)
vars_state(::VorticityModel, ::Prognostic, T) = @vars()
vars_state(::VorticityModel, ::Gradient, T) = @vars()
vars_state(::VorticityModel, ::GradientFlux, T) = @vars()
function init_state_auxiliary!(
    m::VorticityModel,
    state_auxiliary::MPIStateArray,
    grid,
    direction,
) end
function init_state_prognostic!(
    ::VorticityModel,
    state::Vars,
    aux::Vars,
    localgeo,
    t,
) end
function flux_first_order!(
    ::VorticityModel,
    flux::Grad,
    state::Vars,
    auxstate::Vars,
    t::Real,
    direction,
)
    ρ = auxstate.ρ
    ρinv = 1/ρ
    ρu = auxstate.ρu
    u = ρinv * ρu
    @inbounds begin
        flux.ω = @SMatrix [ 0     u[3] -u[2];
                            -u[3]  0     u[1];
                            u[2] -u[1]  0    ]
    end

end
flux_second_order!(::VorticityModel, _...) = nothing
source!(::VorticityModel, _...) = nothing

boundary_conditions(::VorticityModel) = ntuple(i -> nothing, 6)
boundary_state!(nf, ::Nothing, ::VorticityModel, _...) = nothing