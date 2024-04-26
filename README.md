# ManyZonalJets.jl

[![Build Status](https://github.com/SpeedyWeather/ManyZonalJets.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SpeedyWeather/ManyZonalJets.jl/actions/workflows/CI.yml?query=branch%3Amain)

The zonal jet initial conditions from 
[Galewsky, 2004](https://doi.org/10.3402/tellusa.v56i5.14436) but for many/multiple zonal jets for
[SpeedyWeather](https://github.com/SpeedyWeather/SpeedyWeather.jl)'s shallow-water model.

## Usage

```julia
using SpeedyWeather
using ManyZonalJets     # exports ZonalJets

spectral_grid = SpectralGrid(trunc=85, nlev=1)  # T85 resolution, ~165km global

# initial conditions for three jets, at 45˚S, 0˚N, and 45˚N
# all same strength and width, all perturbed the same way
initial_conditions = ZonalJets(latitude=[45,0,-45])
orography = NoOrography(spectral_grid)

# create a model with those initial conditions + initialize
model = ShallowWaterModel(;spectral_grid, initial_conditions)
simulation = initialize!(model)

# run for 6 days
run!(simulation, period=Day(6))
```

## Gallery

## Installation

ManyZonalJets.jl is not registered, so you have to install it manually with
Julia's package manager (open it with `]`)
```julia
(@v1.10) pkg> add https://github.com/SpeedyWeather/ManyZonalJets.jl
```
Or alternatively, `using Pkg`, then `Pkg.add("https://github.com/SpeedyWeather/ManyZonalJets.jl")`
(which is equivalent).