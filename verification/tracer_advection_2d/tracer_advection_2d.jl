using Printf
using Logging
using Plots

using Oceananigans
using Oceananigans.Advection
using Oceananigans.OutputWriters
using Oceananigans.Utils

using Oceananigans.Grids: ynodes, znodes

include("stommel_gyre.jl")

ENV["GKSwstype"] = "100"
pyplot()

Logging.global_logger(OceananigansLogger())

#####
##### Initial (= final) conditions
#####

@inline ϕ_Gaussian(x, y; L, A, σˣ, σʸ) = A * exp(-(x-L/2)^2/(2σˣ^2) -(y-L/2)^2/(2σʸ^2))
@inline ϕ_Square(x, y; L, A, σˣ, σʸ)   = A * (-σˣ <= x-L/2 <= σˣ) * (-σʸ <= y-L/2 <= σʸ)

ic_name(::typeof(ϕ_Gaussian)) = "Gaussian"
ic_name(::typeof(ϕ_Square))   = "Square"

#####
##### Experiment functions
#####

function setup_simulation(N, T, CFL, ϕₐ, advection_scheme; u, v)
    topology = (Flat, Periodic, Periodic)
    domain = (x=(0, 1), y=(0, L), z=(0, L))
    grid = RegularCartesianGrid(topology=topology, size=(1, N, N), halo=(3, 3, 3); domain...)

    model = IncompressibleModel(
             grid = grid,
        advection = advection_scheme,
          tracers = :c,
         buoyancy = nothing,
          closure = IsotropicDiffusivity(ν=0, κ=0)
    )

    set!(model, v=u, w=v, c=ϕₐ)

    v_max = maximum(abs, interior(model.velocities.v))
    w_max = maximum(abs, interior(model.velocities.w))
    Δt = CFL * min(grid.Δy, grid.Δz) / max(v_max, w_max)

    simulation = Simulation(model, Δt=Δt, stop_time=T, progress=print_progress, iteration_interval=10)

    filename = @sprintf("%s_%s_N%d_CFL%.2f.nc", ic_name(ϕₐ), typeof(advection_scheme), N, CFL)
    fields = Dict("v" => model.velocities.v, "w" => model.velocities.w, "c" => model.tracers.c)
    global_attributes = Dict("N" => N, "CFL" => CFL, "advection_scheme" => string(typeof(advection_scheme)))
    output_attributes = Dict("c" => Dict("longname" => "passive tracer"))

    simulation.output_writers[:fields] =
        NetCDFOutputWriter(model, fields, filename=filename, time_interval=0.01,
                           global_attributes=global_attributes, output_attributes=output_attributes)

    return simulation
end

function print_progress(simulation)
    model = simulation.model

    progress = 100 * (model.clock.time / simulation.stop_time)

    v_max = maximum(abs, interior(model.velocities.v))
    w_max = maximum(abs, interior(model.velocities.w))
    c_min, c_max = extrema(interior(model.tracers.c))

    i, t = model.clock.iteration, model.clock.time
    @info @sprintf("[%06.2f%%] i: %d, t: %.4f, U_max: (%.2e, %.2e), c: (min=%.3e, max=%.3e)",
                   progress, i, t, v_max, w_max, c_min, c_max)
end

L = τ₀ = β = 1
r = 0.04β*L

A  = 1
σˣ = σʸ = 0.05L

T = 2

u = (x, y, z) -> 1 # u_Stommel(y, z, L=L, τ₀=τ₀, β=β, r=r)
v = (x, y, z) -> 1 # v_Stommel(y, z, L=L, τ₀=τ₀, β=β, r=r)

ϕ_λ_Gaussian = (x, y, z) -> ϕ_Gaussian(y, z, L=L, A=A, σˣ=σˣ, σʸ=σʸ)
ϕ_λ_Square   = (x, y, z) -> ϕ_Square(y, z, L=L, A=A, σˣ=σˣ, σʸ=σʸ)

ic_name(::typeof(ϕ_λ_Gaussian)) = ic_name(ϕ_Gaussian)
ic_name(::typeof(ϕ_λ_Square))   = ic_name(ϕ_Square)

ϕ_λs = (ϕ_λ_Gaussian, ϕ_λ_Square)
schemes = (CenteredSecondOrder(), CenteredFourthOrder(), WENO5())
Ns = (32, 256)
CFLs = (0.05, 0.3)

for ϕ in ϕ_λs, scheme in schemes, N in Ns, CFL in CFLs
    simulation = setup_simulation(N, T, CFL, ϕ, scheme, u=v, v=v)

    # simulation.stop_time = T/2
    run!(simulation)

    # set!(simulation.model, v=-1, w=-1)

    # simulation.stop_time = T
    # run!(simulation)
end
