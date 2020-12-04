# Insider Model allows computation of variables using model kernels using mini balance laws
using ClimateMachine.BalanceLaws:
    AbstractStateType, Auxiliary, Gradient

import ClimateMachine.BalanceLaws:
    vars_state


# ------------------------ General Insider Model ---------------------- #

vars_state(::InsiderModel, ::AbstractStateType, FT) = @vars()

function flux_first_order!(
    ::InsiderModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    ts,
    direction,
) end
    




# ------------------------ Begin Vorticity Model ---------------------- #
"""
  VorticityModel <: InsiderModel
## References
 - 
"""
struct VorticityModel{FT} <: InsiderModel
    "vertical vorticity component `[1/s]`"
    #ω::FT
end
vars_state(m::VorticityModel, ::Auxiliary, FT) = @vars(vort::FT)

function flux_first_order!(
    m::VorticityModel,
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    ts,
    direction,
)
    ρ = aux.ρ
    ρinv = 1/ρ
    ρu = aux.ρu
    u = ρinv * ρu
    @inbounds begin
    flux.ω = @SMatrix [ 0     u[3] -u[2];
                        -u[3]  0     u[1];
                        u[2] -u[1]  0    ]
    end

end
# -------------------------- End Radiation Model ------------------------ #