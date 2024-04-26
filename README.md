# ManyZonalJets.jl

[![Build Status](https://github.com/SpeedyWeather/ManyZonalJets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SpeedyWeather/ManyZonalJets.jl/actions/workflows/CI.yml?query=branch%3Amain)

The zonal jet initial conditions from 
[Galewsky, 2004](https://doi.org/10.3402/tellusa.v56i5.14436) but for many/multiple zonal jets for
[SpeedyWeather](https://github.com/SpeedyWeather/SpeedyWeather.jl)'s shallow-water model.

## Usage

```julia
using SpeedyWeather
using ManyZonalJets

spectral_grid = SpectralGrid(trunc=85, nlev=1)
initial_conditions = ZonalJets(latitude=[45,0,-45])
orography = NoOrography(spectral_grid)
model = ShallowWaterModel(;spectral_grid, orography, initial_conditions)
simulation = initialize!(model)
run!(simulation, period=Day(6))
```

`using ManyZonalJets` exports the `ZonalJets` type. 
`spectral_grid` defines the resolution to be T85 (~165km global).
`ZonalJets` initial conditions are then created with jets at 45˚N, 0˚, 45˚S.
`NoOrography` is used as in [Galewsky, 2004](https://doi.org/10.3402/tellusa.v56i5.14436).
Then create a `ShallowWaterModel` with those components, `initialize!` and 
`run!` the `simulation` for 6 days.

## Gallery

At T127 vorticity evolves like

https://github.com/SpeedyWeather/ManyZonalJets.jl/assets/25530332/344352f4-ef8a-4017-81d1-8f80da4a5c63

with mountains via `orography = EarthOrography(spectral_grid)` we have instead

https://github.com/SpeedyWeather/ManyZonalJets.jl/assets/25530332/adffaca3-f638-42a0-bd7c-5ca4942535ad

## Installation

ManyZonalJets.jl is not registered, so you have to install it manually with
Julia's package manager (open it with `]`)
```julia
(@v1.10) pkg> add https://github.com/SpeedyWeather/ManyZonalJets.jl
```
Or alternatively, `using Pkg`, then `Pkg.add("https://github.com/SpeedyWeather/ManyZonalJets.jl")`
(which is equivalent).
