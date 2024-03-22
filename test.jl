using Oceananigans

grid = RectilinearGrid(size = (5, 5, 4), 
                       x = (-1e3, 1e3), 
                       y = (-1e3, 1e3), 
                       z = (-1e3, 0),
                       halo = (3,3,3)
                       )

model = HydrostaticFreeSurfaceModel(grid = grid)

ηᵢ(x, y, z) = 1 * exp(-(x^2 + y^2) / 2 * (2e2)^2)

set!(model, η = ηᵢ)

simulation = Simulation(model, Δt=100, stop_time = 1000)

eta = model.free_surface.η

simulation.output_writers[:surface] = JLD2OutputWriter(model, 
                                                       (η=eta,),
                                                       schedule = TimeInterval(200),
                                                       filename = "surface",
                                                       with_halos = false,
                                                       overwrite_existing = true)

# simulation.output_writers[:surface_nc] = NetCDFOutputWriter(model, (η=model.free_surface.η,),
#                                                        schedule = TimeInterval(200),
#                                                        filename = "surface",
#                                                        overwrite_existing = true)