using Oceananigans.Fields: Field, fill_halo_regions!, set!
using Oceananigans.Grids: Grids, bottommost_active_node, AbstractStaticGrid, constructor_arguments
using Oceananigans.Utils: prettysummary

import Oceananigans.Operators: Δrᶜᶜᶜ, Δrᶜᶜᶠ, Δrᶜᶠᶜ, Δrᶜᶠᶠ, Δrᶠᶜᶜ, Δrᶠᶜᶠ, Δrᶠᶠᶜ, Δrᶠᶠᶠ,
                               Δzᶜᶜᶜ, Δzᶜᶜᶠ, Δzᶜᶠᶜ, Δzᶜᶠᶠ, Δzᶠᶜᶜ, Δzᶠᶜᶠ, Δzᶠᶠᶜ, Δzᶠᶠᶠ

#####
##### ZPartialCell
#####

abstract type AbstractGridFitted{H, T} <: AbstractGridFittedBoundary end

struct ZPartialCell{H, T, E} <: AbstractGridFitted{H, T}
    bottom_height :: H
    top_height :: T
    minimum_fractional_cell_height :: E
end

const PartialCellBottom{T} = ZPartialCell{<:Any, Nothing, T} where T 

const PCBIBG{FT, TX, TY, TZ} = ImmersedBoundaryGrid{FT, TX, TY, TZ, <:Any, <:ZPartialCell} where {FT, TX, TY, TZ}

function Base.summary(ib::ZPartialCell)
    bzmax = maximum(parent(ib.bottom_height))
    bzmin = minimum(parent(ib.bottom_height))
    bzmean = mean(parent(ib.bottom_height))

    
    summary1 = "PartialCellBottom \n"

    summary2 = string("│   ├── bottom_height => ",
                      "mean(zb)=", prettysummary(bzmean),
                      ", min(zb)=", prettysummary(bzmin),
                      ", max(zb)=", prettysummary(bzmax),
                      "\n")


    if isnothing(ib.top_height)
        summary3 = string("│   ├── top_height => ", "None", "\n")
    else
        summary3 = string(
            "│   ├── top_height => ",
            "mean(zt)=$(prettysummary(mean(parent(ib.top_height)))),",
            "min(zt)=$(prettysummary(minimum(parent(ib.top_height)))),", 
            "max(zt)=$(prettysummary(maximum(parent(ib.top_height))))"
        )
    end
    
    summary4 = string("│   └── minimum_fractional_cell_height => ",
                      "ϵ=$(prettysummary(ib.minimum_fractional_cell_height))")
        
    return summary1 * summary2 * summary3 * summary4
end

Base.summary(ib::ZPartialCell{<:Function}) = @sprintf("ZPartialCell(%s, ϵ=%.1f)",
                                                           prettysummary(ib.bottom_height, false),
                                                           ib.minimum_fractional_cell_height)

function Base.show(io::IO, ib::ZPartialCell)
    print(io, summary(ib), '\n')
    print(io, "├── bottom_height: ", prettysummary(ib.bottom_height), '\n')
    print(io, "├── top_height: ", prettysummary(ib.top_height), '\n')
    print(io, "└── minimum_fractional_cell_height: ", prettysummary(ib.minimum_fractional_cell_height))
end

"""
    PartialCellBottom(bottom_height; minimum_fractional_cell_height=0.2)

Return `PartialCellBottom` representing an immersed boundary with "partial"
bottom cells. That is, the height of the bottommost cell in each column is reduced
to fit the provided `bottom_height`, which may be a `Field`, `Array`, or function
of `(x, y)`.

The height of partial bottom cells is greater than

```
minimum_fractional_cell_height * Δz,
```

where `Δz` is the original height of the bottom cell underlying grid.
"""
function ZPartialCell(bottom_height, top_height; minimum_fractional_cell_height=0.2)
    return ZPartialCell(bottom_height, top_height, minimum_fractional_cell_height)
end

function ZPartialCell(bottom_height; minimum_fractional_cell_height=0.2)
    return ZPartialCell(bottom_height, nothing, minimum_fractional_cell_height)
end

function materialize_immersed_boundary(grid, ib::ZPartialCell)
    bottom_field = Field{Center, Center, Nothing}(grid)
    set!(bottom_field, ib.bottom_height)

    if !isnothing(ib.top_height)
        top_field = Field{Center, Center, Nothing}(grid)
        set!(top_field, ib.top_height)
    else
        top_field = nothing
    end

    minimum_fractional_cell_height = convert(eltype(grid), ib.minimum_fractional_cell_height)
    new_ib = ZPartialCell(bottom_field, top_field, minimum_fractional_cell_height)

    @apply_regionally compute_numerical_bottom_height!(bottom_field, grid, new_ib)
    fill_halo_regions!(bottom_field)

    return new_ib
end

@kernel function _compute_numerical_bottom_height!(bottom_field, grid, ib::ZPartialCell)
    i, j = @index(Global, NTuple)

    # Save analytical bottom height
    zb = @inbounds bottom_field[i, j, 1]

    # Cap bottom height at Lz and at rnode(i, j, grid.Nz+1, grid, c, c, f)

    domain_bottom = rnode(i, j, 1, grid, c, c, f)
    domain_top    = rnode(i, j, grid.Nz+1, grid, c, c, f)
    @inbounds bottom_field[i, j, 1] = clamp(zb, domain_bottom, domain_top)
    adjusted_zb = bottom_field[i, j, 1]

    ϵ  = ib.minimum_fractional_cell_height

    for k in 1:grid.Nz
        z⁻ = rnode(i, j, k,   grid, c, c, f)
        z⁺ = rnode(i, j, k+1, grid, c, c, f)
        Δz = Δrᶜᶜᶜ(i, j, k, grid)
        bottom_cell = z⁻ ≤ adjusted_zb < z⁺
        capped_zb   = min(z⁺ - ϵ * Δz, adjusted_zb)

        # If the size of the bottom cell is less than ϵ Δz,
        # we enforce a minimum size of ϵ Δz.
        adjusted_zb = ifelse(bottom_cell, capped_zb, adjusted_zb)
    end
    @inbounds bottom_field[i, j, 1] = adjusted_zb
end

function Architectures.on_architecture(arch, ib::ZPartialCell{<:Field})
    architecture(ib.bottom_height) == arch && return ib
    arch_grid = on_architecture(arch, ib.bottom_height.grid)
    new_bottom_height = Field{Center, Center, Nothing}(arch_grid)
    copyto!(parent(new_bottom_height), parent(ib.bottom_height))
    return ZPartialCell(new_bottom_height, ib.minimum_fractional_cell_height)
end

Adapt.adapt_structure(to, ib::ZPartialCell) = ZPartialCell(adapt(to, ib.bottom_height),
                                                                     ib.minimum_fractional_cell_height)

Architectures.on_architecture(to, ib::ZPartialCell) = ZPartialCell(on_architecture(to, ib.bottom_height),
                                                                             on_architecture(to, ib.minimum_fractional_cell_height))

"""
    immersed     underlying

      --x--        --x--


        ∘   ↑        ∘   k+1
            |
            |
  k+1 --x-- |  k+1 --x--    ↑      <- node z
        ∘   ↓               |
   zb ⋅⋅x⋅⋅                 |
                            |
                     ∘   k  | Δz
                            |
                            |
                 k --x--    ↓

Criterion is zb ≥ z - ϵ Δz

"""
@inline function _immersed_cell(i, j, k, underlying_grid, ib::ZPartialCell)
    r⁺ = rnode(i, j, k + 1, underlying_grid, c, c, f)
    ϵ  = ib.minimum_fractional_cell_height
    Δr = Δrᶜᶜᶜ(i, j, k, underlying_grid)
    r★ = r⁺ - Δr * ϵ
    rᵇ = @inbounds ib.bottom_height[i, j, 1]
    return r★ < rᵇ
end

@inline function Δrᶜᶜᶜ(i, j, k, ibg::PCBIBG)
    underlying_grid = ibg.underlying_grid
    ib = ibg.immersed_boundary

    # Get node at face above and defining nodes on c,c,f
    r⁺ = rnode(i, j, k + 1, underlying_grid, c, c, f)

    # Get bottom r-coordinate and fractional Δr parameter
    rᵇ = @inbounds ib.bottom_height[i, j, 1]

    # Are we in a bottom cell?
    at_the_bottom = bottommost_active_node(i, j, k, ibg, c, c, c)

    full_Δr    = Δrᶜᶜᶜ(i, j, k, ibg.underlying_grid)
    partial_Δr = r⁺ - rᵇ

    return ifelse(at_the_bottom, partial_Δr, full_Δr)
end

@inline function Δrᶜᶜᶠ(i, j, k, ibg::PCBIBG)
    just_above_bottom = bottommost_active_node(i, j, k, ibg, c, c, f)
    rᶜ = rnode(i, j, k, ibg.underlying_grid, c, c, c)
    rᶠ = rnode(i, j, k, ibg.underlying_grid, c, c, f)

    full_Δr    = Δrᶜᶜᶠ(i, j, k, ibg.underlying_grid)
    partial_Δr = rᶜ - rᶠ + Δrᶜᶜᶜ(i, j, k-1, ibg) / 2

    return ifelse(just_above_bottom, partial_Δr, full_Δr)
end

@inline Δrᶠᶜᶜ(i, j, k, ibg::PCBIBG) = min(Δrᶜᶜᶜ(i-1, j, k, ibg), Δrᶜᶜᶜ(i, j, k, ibg))
@inline Δrᶜᶠᶜ(i, j, k, ibg::PCBIBG) = min(Δrᶜᶜᶜ(i, j-1, k, ibg), Δrᶜᶜᶜ(i, j, k, ibg))
@inline Δrᶠᶠᶜ(i, j, k, ibg::PCBIBG) = min(Δrᶠᶜᶜ(i, j-1, k, ibg), Δrᶠᶜᶜ(i, j, k, ibg))

@inline Δrᶠᶜᶠ(i, j, k, ibg::PCBIBG) = min(Δrᶜᶜᶠ(i-1, j, k, ibg), Δrᶜᶜᶠ(i, j, k, ibg))
@inline Δrᶜᶠᶠ(i, j, k, ibg::PCBIBG) = min(Δrᶜᶜᶠ(i, j-1, k, ibg), Δrᶜᶜᶠ(i, j, k, ibg))
@inline Δrᶠᶠᶠ(i, j, k, ibg::PCBIBG) = min(Δrᶠᶜᶠ(i, j-1, k, ibg), Δrᶠᶜᶠ(i, j, k, ibg))

# Make sure Δz works for horizontally-Flat topologies.
# (There's no point in using z-Flat with ZPartialCell).
XFlatPCBIBG = ImmersedBoundaryGrid{<:Any, <:Flat, <:Any, <:Any, <:Any, <:ZPartialCell}
YFlatPCBIBG = ImmersedBoundaryGrid{<:Any, <:Any, <:Flat, <:Any, <:Any, <:ZPartialCell}

@inline Δrᶠᶜᶜ(i, j, k, ibg::XFlatPCBIBG) = Δrᶜᶜᶜ(i, j, k, ibg)
@inline Δrᶠᶜᶠ(i, j, k, ibg::XFlatPCBIBG) = Δrᶜᶜᶠ(i, j, k, ibg)
@inline Δrᶜᶠᶜ(i, j, k, ibg::YFlatPCBIBG) = Δrᶜᶜᶜ(i, j, k, ibg)

@inline Δrᶜᶠᶠ(i, j, k, ibg::YFlatPCBIBG) = Δrᶜᶜᶠ(i, j, k, ibg)
@inline Δrᶠᶠᶜ(i, j, k, ibg::XFlatPCBIBG) = Δrᶜᶠᶜ(i, j, k, ibg)
@inline Δrᶠᶠᶜ(i, j, k, ibg::YFlatPCBIBG) = Δrᶠᶜᶜ(i, j, k, ibg)

# Vertically-static, partial cell bottom, immersed boundary grid
VSPCBIBG = ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:AbstractStaticGrid, <:ZPartialCell}
@inline Δzᶜᶜᶜ(i, j, k, ibg::VSPCBIBG) = Δrᶜᶜᶜ(i, j, k, ibg)
@inline Δzᶠᶜᶜ(i, j, k, ibg::VSPCBIBG) = Δrᶠᶜᶜ(i, j, k, ibg)
@inline Δzᶜᶠᶜ(i, j, k, ibg::VSPCBIBG) = Δrᶜᶠᶜ(i, j, k, ibg)
@inline Δzᶜᶜᶠ(i, j, k, ibg::VSPCBIBG) = Δrᶜᶜᶠ(i, j, k, ibg)
@inline Δzᶠᶠᶜ(i, j, k, ibg::VSPCBIBG) = Δrᶠᶠᶜ(i, j, k, ibg)
@inline Δzᶜᶠᶠ(i, j, k, ibg::VSPCBIBG) = Δrᶜᶠᶠ(i, j, k, ibg)
@inline Δzᶠᶜᶠ(i, j, k, ibg::VSPCBIBG) = Δrᶠᶜᶠ(i, j, k, ibg)
@inline Δzᶠᶠᶠ(i, j, k, ibg::VSPCBIBG) = Δrᶠᶠᶠ(i, j, k, ibg)

function Grids.constructor_arguments(grid::PCBIBG)
    underlying_grid_args, underlying_grid_kwargs = constructor_arguments(grid.underlying_grid)
    partial_cell_bottom_args = Dict(:bottom_height => grid.immersed_boundary.bottom_height,
                                    :minimum_fractional_cell_height => grid.immersed_boundary.minimum_fractional_cell_height)
    return underlying_grid_args, underlying_grid_kwargs, partial_cell_bottom_args
end

function Base.:(==)(pcb1::ZPartialCell, pcb2::ZPartialCell)
    return pcb1.bottom_height == pcb2.bottom_height && pcb1.minimum_fractional_cell_height == pcb2.minimum_fractional_cell_height
end
