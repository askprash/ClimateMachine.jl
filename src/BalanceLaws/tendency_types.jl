#### Tendency types

# Terminology:
#
# `∂_t Yᵢ + (∇•F₁(Y))ᵢ + (∇•F₂(Y,G)))ᵢ = (S(Y,G))ᵢ`
# `__Tᵤ__   ____T₁____   ______T₂______   ___S___`
#
#  - `Yᵢ` - the i-th prognostic variable
#  - `Y` - the prognostic state (column vector)
#  - `G = ∇Y` - the gradient of the prognostic state (rank 2 tensor)
#  - `F₁` - the first order tendency flux (rank 2 tensor)
#  - `F₂` - the second order tendency flux (rank 2 tensor)
#  - `Tᵤ` - the unsteady tendency (column vector)
#  - `T₁` - the first order flux tendency (column vector)
#  - `T₂` - the second order flux tendency (column vector)
#  - `S` - the non-conservative tendency source (column vector)

export PrognosticVariable
export FirstOrder, SecondOrder
export Flux, Source
export TendencyDef

"""
    PrognosticVariable

Subtypes are used for specifying
each prognostic variable.
"""
abstract type PrognosticVariable end

"""
    AbstractOrder

Subtypes are used for dispatching
on the flux order.
"""
abstract type AbstractOrder end
struct FirstOrder <: AbstractOrder end
struct SecondOrder <: AbstractOrder end

"""
    AbstractTendencyType

Subtypes are used for specifying a
tuple of tendencies to be accumulated.
"""
abstract type AbstractTendencyType end
struct Flux{O <: AbstractOrder} <: AbstractTendencyType end
struct Source <: AbstractTendencyType end

"""
    TendencyDef

Subtypes are used for specifying
each tendency definition.
"""
abstract type TendencyDef{TT <: AbstractTendencyType, PV <: PrognosticVariable} end

"""
    eq_tends(::BalanceLaw, ::AbstractTendencyType)

Calls `eq_tends(::Tuple, ::BalanceLaw, ::AbstractTendencyType)`.
All of the tendencies for the given balance law

    eq_tends(::Tuple, ::BalanceLaw, ::AbstractTendencyType)

Calls `eq_tends(::PrognosticVariable, ::BalanceLaw, ::AbstractTendencyType)`

    eq_tends(::PrognosticVariable, ::BalanceLaw, ::AbstractTendencyType)

A tuple of `TendencyDef`s given
 - `PrognosticVariable` prognostic variable
 - `AbstractTendencyType` tendency type
 - `BalanceLaw` balance law

i.e., a tuple of `TendencyDef`s corresponding
to `F₁`, `F₂`, **or** `S` for a single
prognostic variable in:

    `∂_t Yᵢ + (∇•F₁(Y))ᵢ + (∇•F₂(Y,G)))ᵢ = (S(Y,G))ᵢ`
"""
function eq_tends end

function eq_tends(bl::BalanceLaw, tt::AbstractTendencyType)
    tup = ()
    map(prognostic_vars(bl)) do pv
        eqt = eq_tends(pv, bl, tt)
        if !isempty(eqt)
            tup = (tup..., eqt...)
        end
    end
    for p in propertynames(bl)
        subbl = getproperty(bl, p)
        subbl_pvs = prognostic_vars(subbl)
        isempty(subbl_pvs) && continue
        subeqt = map(subbl_pvs) do pv
            eqt = eq_tends(pv, subbl, tt)
            if !isempty(eqt)
                tup = (tup..., eqt...)
            end
        end
    end
    return tup
end

"""
    prognostic_vars(::Any)

A tuple of `PrognosticVariable`s given
the `BalanceLaw`.

i.e., a tuple of `PrognosticVariable`s
corresponding to a _part_ of the column-
vector `Yᵢ` in:

    `∂_t Yᵢ + (∇•F₁(Y))ᵢ + (∇•F₂(Y,G)))ᵢ = (S(Y,G))ᵢ`
"""
prognostic_vars(::Any) = ()

# Skip functions/reals AbstractArrays
all_prognostic_vars(bl::T, tup) where {T<:Function} = tup
all_prognostic_vars(bl::T, tup) where {T<:Real} = tup
all_prognostic_vars(bl::T, tup) where {T<:AbstractArray} = tup

"""
    all_prognostic_vars(::BalanceLaw)

A tuple of `PrognosticVariable`s given
the `BalanceLaw`.

i.e., a tuple of `PrognosticVariable`s
corresponding to the _entire_ column-
vector `Yᵢ` in:

    `∂_t Yᵢ + (∇•F₁(Y))ᵢ + (∇•F₂(Y,G)))ᵢ = (S(Y,G))ᵢ`

by default, we just call `prognostic_vars`
"""
function all_prognostic_vars(bl, tup = ())
    pvs = prognostic_vars(bl)
    if !isempty(pvs)
        tup = (tup..., pvs...)
    end
    pns = propertynames(bl)
    if !isempty(pns)
        for p in propertynames(bl)
            tup = all_prognostic_vars(getproperty(bl, p), tup)
        end
    end
    return tup
end

export sources
"""
    sources(bl::BalanceLaw)

A tuple of `TendencyDef{Source}`s
given the `BalanceLaw`.

i.e., a tuple of `TendencyDef{Source}`s
corresponding to the column-vector `S` in:

    `∂_t Yᵢ + (∇•F₁(Y))ᵢ + (∇•F₂(Y,G)))ᵢ = (S(Y,G))ᵢ`
"""
function sources(bl::BalanceLaw)
    tend = eq_tends.(all_prognostic_vars(bl), Ref(bl), Ref(Source()))
    tend = filter(x -> x ≠ nothing, tend)
    return Tuple(Iterators.flatten(tend))
end

export fluxes
"""
    fluxes(bl::BalanceLaw, order::O) where {O <: AbstractOrder}

A tuple of `TendencyDef{Flux{O}}`s
given the `BalanceLaw` and the `order::O`.

i.e., a tuple of `TendencyDef{Flux{O}}`s
corresponding to the column-vector `F₁`
or `F₂` given the flux order `order::O` in:

    `∂_t Yᵢ + (∇•F₁(Y))ᵢ + (∇•F₂(Y,G)))ᵢ = (S(Y,G))ᵢ`
"""
function fluxes(bl::BalanceLaw, order::O) where {O <: AbstractOrder}
    tend = eq_tends.(all_prognostic_vars(bl), Ref(bl), Ref(Flux{O}()))
    tend = filter(x -> x ≠ nothing, tend)
    return Tuple(Iterators.flatten(tend))
end
