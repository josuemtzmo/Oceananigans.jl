using Oceananigans.Fields: validate_indices, Reduction
using Oceananigans.AbstractOperations: AbstractOperation, ComputedField
using Oceananigans.Grids: default_indices

restrict_to_interior(::Colon, loc, topo, N, H) = interior_indices(loc, topo, N)
restrict_to_interior(::Colon, ::Nothing, topo, N, H) = UnitRange(1, 1)
restrict_to_interior(index::UnitRange, ::Nothing, topo, N, H) = UnitRange(1, 1)

function restrict_to_interior(index::UnitRange, loc, topo, N, H)
    from = max(first(index), 1)
    to = min(last(index), last(interior_indices(loc, topo, N)))
    if from == to
        from = 1 - H
        to   = 1 - H
    end
    return UnitRange(from, to)
end

#####
##### Function output fallback
#####

function construct_output(output, grid, indices, with_halos)
    if !(indices isa typeof(default_indices(3)))
        output_type = output isa Function ? "Function" : ""
        @warn "Cannot slice $output_type $output with $indices: output will be unsliced."
    end

    return output
end

#####
##### Support for Field, Reduction, and AbstractOperation outputs
#####

function output_indices(output::Union{AbstractField, Reduction}, grid, indices, with_halos)
    indices = validate_indices(output.indices, location(output), grid)

    if !with_halos # Maybe chop those indices
        loc = map(instantiate, location(output))
        topo = map(instantiate, topology(grid))
        H = (grid.Hx,grid.Hy,grid.Hz)
        indices = map(restrict_to_interior, indices, loc, topo, size(grid),H)
    end

    return indices
end

function construct_output(user_output::Union{AbstractField, Reduction}, grid, user_indices, with_halos)
    indices = output_indices(user_output, grid, user_indices, with_halos)
    return Field(user_output; indices)
end

#####
##### Time-averaging
#####

function construct_output(averaged_output::WindowedTimeAverage{<:Field}, grid, indices, with_halos)
    output = construct_output(averaged_output.operand, grid, indices, with_halos)
    return WindowedTimeAverage(output; schedule=averaged_output.schedule)
end

