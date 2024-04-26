using SpeedyWeather, ManyZonalJets
using Test

@testset "Various jets" begin
    spectral_grid = SpectralGrid(trunc=31, nlev=1)
    initial_conditions = ZonalJets()
    initial_conditions = ZonalJets(latitude=[45,0,-45])
    initial_conditions = ZonalJets(latitude=[45,0,-45], umax=[80, 40, 20])
    initial_conditions = ZonalJets(latitude=[45,0,-45], perturb_height=[120, 0, 0])
    model = ShallowWaterModel(;spectral_grid, initial_conditions)
    simulation = initialize!(model)
    run!(simulation, period=Day(6))
end
