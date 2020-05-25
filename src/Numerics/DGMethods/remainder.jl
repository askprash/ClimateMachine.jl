using StaticNumbers
export remainder_DGModel

import ..BalanceLaws:
    vars_state_conservative,
    vars_state_auxiliary,
    vars_state_gradient,
    vars_state_gradient_flux,
    vars_integrals,
    vars_reverse_integrals,
    vars_gradient_laplacian,
    vars_hyperdiffusive,
    init_state_conservative!,
    init_state_auxiliary!,
    flux_first_order!,
    flux_second_order!,
    compute_gradient_flux!,
    compute_gradient_argument!,
    source!,
    transform_post_gradient_laplacian!,
    wavespeed,
    boundary_state!,
    update_auxiliary_state!,
    integral_load_auxiliary_state!,
    integral_set_auxiliary_state!,
    reverse_integral_load_auxiliary_state!,
    reverse_integral_set_auxiliary_state!

"""
    RemainderModel(main::BalanceLaw, subcomponents::Tuple)

Compute the "remainder" contribution of the `main` model, after subtracting
`subcomponents`.

Currently only the `flux_nondiffusive!` and `source!` are handled by the
remainder model
"""
struct RemainderModel{M, S} <: BalanceLaw
    main::M
    subs::S
end

function remainder_DGModel(
    maindg,
    subsdg,
    numerical_flux_first_order,
    numerical_flux_second_order,
    numerical_flux_gradient;
    state_auxiliary = maindg.state_auxiliary,
    state_gradient_flux = create_gradient_state(
        maindg.balance_law,
        maindg.grid,
    ),
    states_higher_order = create_higher_order_states(
        maindg.balance_law,
        maindg.grid,
    ),
    direction = EveryDirection(),
    diffusion_direction = maindg.diffusion_direction,
    modeldata = maindg.modeldata,
)
    balance_law = RemainderModel(
        maindg.balance_law,
        ntuple(i -> subsdg[i].balance_law, length(subsdg)),
    )
    DGModel(
        balance_law,
        maindg.grid,
        numerical_flux_first_order,
        numerical_flux_second_order,
        numerical_flux_gradient,
        state_auxiliary,
        state_gradient_flux,
        states_higher_order,
        direction,
        diffusion_direction,
        modeldata,
    )
end

# Inherit most of the functionality from the main model
vars_state_conservative(rem_balance_law::RemainderModel, FT) =
    vars_state_conservative(rem_balance_law.main, FT)

vars_state_gradient(rem_balance_law::RemainderModel, FT) =
    vars_state_gradient(rem_balance_law.main, FT)

vars_state_gradient_flux(rem_balance_law::RemainderModel, FT) =
    vars_state_gradient_flux(rem_balance_law.main, FT)

vars_state_auxiliary(rem_balance_law::RemainderModel, FT) =
    vars_state_auxiliary(rem_balance_law.main, FT)

vars_integrals(rem_balance_law::RemainderModel, FT) =
    vars_integrals(rem_balance_law.main, FT)

vars_reverse_integrals(rem_balance_law::RemainderModel, FT) =
    vars_integrals(rem_balance_law.main, FT)

vars_gradient_laplacian(rem_balance_law::RemainderModel, FT) =
    vars_gradient_laplacian(rem_balance_law.main, FT)

vars_hyperdiffusive(rem_balance_law::RemainderModel, FT) =
    vars_hyperdiffusive(rem_balance_law.main, FT)

update_auxiliary_state!(dg::DGModel, rem_balance_law::RemainderModel, args...) =
    update_auxiliary_state!(dg, rem_balance_law.main, args...)

update_auxiliary_state_gradient!(
    dg::DGModel,
    rem_balance_law::RemainderModel,
    args...,
) = update_auxiliary_state_gradient!(dg, rem_balance_law.main, args...)

integral_load_auxiliary_state!(rem_balance_law::RemainderModel, args...) =
    integral_load_auxiliary_state!(rem_balance_law.main, args...)

integral_set_auxiliary_state!(rem_balance_law::RemainderModel, args...) =
    integral_set_auxiliary_state!(rem_balance_law.main, args...)

reverse_integral_load_auxiliary_state!(
    rem_balance_law::RemainderModel,
    args...,
) = reverse_integral_load_auxiliary_state!(rem_balance_law.main, args...)

reverse_integral_set_auxiliary_state!(
    rem_balance_law::RemainderModel,
    args...,
) = reverse_integral_set_auxiliary_state!(rem_balance_law.main, args...)

transform_post_gradient_laplacian!(rem_balance_law::RemainderModel, args...) =
    transform_post_gradient_laplacian!(rem_balance_law.main, args...)

flux_second_order!(rem_balance_law::RemainderModel, args...) =
    flux_second_order!(rem_balance_law.main, args...)

compute_gradient_argument!(rem_balance_law::RemainderModel, args...) =
    compute_gradient_argument!(rem_balance_law.main, args...)

compute_gradient_flux!(rem_balance_law::RemainderModel, args...) =
    compute_gradient_flux!(rem_balance_law.main, args...)

boundary_state!(nf, rem_balance_law::RemainderModel, args...) =
    boundary_state!(nf, rem_balance_law.main, args...)

init_state_auxiliary!(rem_balance_law::RemainderModel, args...) =
    init_state_auxiliary!(rem_balance_law.main, args...)

init_state_conservative!(rem_balance_law::RemainderModel, args...) =
    init_state_conservative!(rem_balance_law.main, args...)

function wavespeed(rem::RemainderModel, nM, state::Vars, aux::Vars, t::Real)
    FT = eltype(state)

    ws = fill(-zero(FT), MVector{number_state_conservative(rem.main, FT), FT})
    rs = fill(-zero(FT), MVector{number_state_conservative(rem.main, FT), FT})

    ws .= wavespeed(rem.main, nM, state, aux, t)

    for sub in rem.subs
        num_state = static(number_state_conservative(sub, Float32))
        @inbounds rs[static(1):num_state] .+= wavespeed(sub, nM, state, aux, t)
    end

    ws .-= rs

    return ws
end

import .NumericalFluxes: normal_boundary_flux_second_order!
boundary_state!(nf, rem_balance_law::RemainderModel, x...) =
    boundary_state!(nf, rem_balance_law.main, x...)

normal_boundary_flux_second_order!(
    nf,
    rem_balance_law::RemainderModel,
    fluxᵀn::Vars{S},
    args...,
) where {S} = normal_boundary_flux_second_order!(
    nf,
    rem_balance_law.main,
    fluxᵀn,
    args...,
)

init_state_auxiliary!(rem_balance_law::RemainderModel, _...) = nothing
init_state_conservative!(rem_balance_law::RemainderModel, _...) = nothing

function flux_first_order!(
    rem_balance_law::RemainderModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    m = parent(flux)
    flux_first_order!(rem_balance_law.main, flux, state, aux, t)

    flux_s = similar(flux)
    m_s = parent(flux_s)

    for sub in rem_balance_law.subs
        fill!(m_s, -zero(eltype(m_s)))
        flux_first_order!(sub, flux_s, state, aux, t)
        m .-= m_s
    end
    nothing
end

function source!(
    rem_balance_law::RemainderModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    m = parent(source)
    source!(rem_balance_law.main, source, state, diffusive, aux, t, direction)

    source_s = similar(source)
    m_s = parent(source_s)

    for sub in rem_balance_law.subs
        fill!(m_s, -zero(eltype(m_s)))
        source!(sub, source_s, state, diffusive, aux, t, direction)
        m .-= m_s
    end
    nothing
end
