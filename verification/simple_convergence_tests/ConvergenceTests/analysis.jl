function extract_two_solutions(analytical_solution, filename)
    grid = RegularCartesianGrid(filename)
    iters = iterations(filename)

    @show filename grid

    u_data = field_data(filename, :u, iters[end])
    u_simulation = Field{Face, Cell, Cell}(u_data, grid, 
                                           FieldBoundaryConditions(grid, (Face, Cell, Cell)))

    Nx, Ny, Nz = size(grid)

    # Works on all grids.
    u_simulation = interior(u_simulation)[1:Nx, 1:Ny, 1:Nz]

    XU = repeat(grid.xF[1:Nx], 1, Ny, Nz)
    YU = repeat(reshape(grid.yC, 1, Ny, 1), Nx, 1, Nz)
    ZU = repeat(reshape(grid.zC, 1, 1, Nz), Nx, Ny, 1)

    t = iteration_time(filename, iters[end])

    u_analytical = analytical_solution.(XU, YU, ZU, t)

    return u_simulation, u_analytical
end

function compute_error(u_simulation, u_analytical)
    absolute_error = @. abs(u_simulation - u_analytical)
    L₁ = mean(absolute_error)
    L∞ = maximum(absolute_error)

    return (L₁=L₁, L∞=L∞)
end

function compute_error(analytical_solution, filename::String)
    u_simulation, u_analytical = extract_two_solutions(analytical_solution, filename)
    return compute_error(u_simulation, u_analytical)
end

compute_errors(analytical_solution, filenames::String...) = 
    [compute_error(analytical_solution, filename) for filename in filenames]

extract_sizes(filenames...) = [size(RegularCartesianGrid(filename)) for filename in filenames]