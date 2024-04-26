module ManyZonalJets

using SpeedyWeather
using DocStringExtensions

export ZonalJets

"""A struct that contains all parameters for multiple zonal jets each
following Galewsky et al, 2004 for intitial conditions for the shallow water model.
Default values as in Galewsky.

All inputs need to be vectors of the same length, but can be omitted
in which case the defaults are applied to each element = each jet.
Only `latitude` always has to be provided for >1 jets as its length
is used to determine how many jets there will be.

    ZonalJets(latitude=[45,0,-45])      # three jets at 45˚N, 0˚, 45˚S
    ZonalJets(  latitude=[45,0,-45],    
                umax=[80, 40, 20])      # makes the southern jets weaker [m/s]
    
The initial conditions for u and η for each jet are superimposed.
E.g. two jets at 45˚N at 80m/s will yield one jet at 160m/s.
Which can be an unstable initial condition. Fields are

$(TYPEDFIELDS)"""
@kwdef mutable struct ZonalJets <: SpeedyWeather.AbstractInitialConditions
    "[OPTION] jet latitude [˚N]"
    latitude::Vector{Float64} = [45]

    "Number of jets, default obtained from length(latitude)"
    njets::Int = length(latitude)
    
    "[OPTION] jet width [˚], default ≈ 19.29˚"
    width::Vector{Float64} = fill((1/4-1/7)*180,njets)

    "[OPTION] jet maximum velocity [m/s]"
    umax::Vector{Float64} = fill(80, njets)
    
    "[OPTION] perturbation latitude [˚N], position in jet by default"
    perturb_lat::Vector{Float64} = latitude
    
    "[OPTION] perturbation longitude [˚E]"
    perturb_lon::Vector{Float64} = fill(270, njets)
    
    "[OPTION] perturbation zonal extent [˚], default ≈ 19.1˚"
    perturb_xwidth::Vector{Float64} = fill(1/3*360/2π, njets)

    "[OPTION] perturbation meridinoal extent [˚], default ≈ 3.8˚"
    perturb_ywidth::Vector{Float64} = fill(1/15*360/2π, njets)
    
    "[OPTION] perturbation amplitude [m]"
    perturb_height::Vector{Float64} = fill(120, njets)
end

"""
$(TYPEDSIGNATURES)
Initial conditions from Galewsky, 2004, Tellus but for each jet in
`ManyJets` and each perturbation in the interface height η superimposed.
Solves for balanced initial conditions (in the case of no perturbation)
for all jets."""
function SpeedyWeather.initialize!(   
    progn::PrognosticVariables,
    initial_conditions::ZonalJets,
    model::ShallowWater
)

    (; latitude, width, umax) = initial_conditions               # for jet
    (; perturb_lat, perturb_lon, perturb_xwidth,                 # for perturbation
        perturb_ywidth, perturb_height) = initial_conditions

    θ₀s = @. (latitude-width)/360*2π    # southern boundary of jet [radians]
    θ₁s = @. (latitude+width)/360*2π    # northern boundary of jet
    eₙs = @. exp(-4/(θ₁s-θ₀s)^2)        # normalisation
    
    θ₂s = perturb_lat*2π/360            # perturbation latitude [radians]
    αs = perturb_xwidth*2π/360          # zonal extent of interface perturbation [radians]
    βs = perturb_ywidth*2π/360          # meridional extent of interface perturbation [radians]
    λs = perturb_lon*2π/360             # perturbation longitude [radians]

    (; rotation, gravity) = model.planet
    (; Grid, NF, radius, nlat_half) = model.spectral_grid
    (; coslat⁻¹) = model.geometry

    u_grid = zeros(Grid{NF}, nlat_half)
    η_perturb_grid = zeros(Grid{NF}, nlat_half)
    lat = RingGrids.get_lat(Grid, nlat_half)
    _, lons = RingGrids.get_colatlons(Grid, nlat_half)

    # loop over every jet and accumulate zonal velocities u and perturbations for η
    for (umaxjet, θ₀, θ₁, eₙ, θ₂, α, β, λ, p) in zip(umax, θ₀s, θ₁s, eₙs, θ₂s, αs, βs, λs, perturb_height)
        for (j, ring) in enumerate(eachring(u_grid, η_perturb_grid))
            θ = lat[j]             # latitude in radians
            
            # velocity per latitude
            if θ₀ < θ < θ₁
                u_θ = umaxjet/eₙ*exp(1/(θ-θ₀)/(θ-θ₁))  # u as in Galewsky, 2004
            else
                u_θ = 0
            end

            # lon-constant part of perturbation
            ηθ = p*cos(θ)*exp(-((θ₂-θ)/β)^2)

            # store in all longitudes
            for ij in ring
                # accumulate here to superimpose all jets
                u_grid[ij] += u_θ/radius*coslat⁻¹[j]   # include scaling for curl!
                
                # calculate perturbation (possibly shifted in lon compared to Galewsky 2004)
                ϕ = lons[ij] - λ
                # accumulate here to superimpose all perturbations
                η_perturb_grid[ij] += exp(-(ϕ/α)^2)*ηθ
            end
        end
    end

    # the following obtain initial conditions for η from u, v=0 via
    # 0 = -∇⋅((ζ+f)*(-v, u)) - ∇²((u^2 + v^2)/2), i.e.
    # invert the Laplacian for
    # ∇²(gη) = -∇⋅(0, (ζ+f)*u) - ∇²(u^2/2)
    u = spectral(u_grid, model.spectral_transform)
    
    # get vorticity initial conditions from curl of u, v
    v = zero(u)     # meridional velocity zero for these initial conditions
    (; vor) = progn.layers[end].timesteps[1]
    curl!(vor, u, v, model.spectral_transform)

    # compute the div = -∇⋅(0,(ζ+f)*u) = -∇×((ζ+f)*u, 0) term, v=0
    vor_grid = gridded(vor, model.spectral_transform)
    f = coriolis(vor_grid; rotation)

    # includes 1/coslat/radius from above for curl!
    # but *radius^2 for the ∇⁻²! operation below!
    vor_flux_grid = @. (vor_grid + f) * u_grid * radius^2
    vor_flux = spectral(vor_flux_grid, model.spectral_transform)
    div = zero(v)
    curl!(div, vor_flux, v, model.spectral_transform)

    # compute the -∇²(u^2/2) term, add to div, divide by gravity
    RingGrids.scale_coslat!(u_grid)     # remove coslat scaling
    u_grid .*= radius                   # no radius scaling as we'll apply ∇⁻²(∇²) (would cancel)
    @. u_grid = convert(NF,1/2) * u_grid^2
    u²_half = SpeedyTransforms.spectral!(u, u_grid, model.spectral_transform)
    ∇²!(div, u²_half, model.spectral_transform, flipsign=true, add=true)
    div .*= inv(gravity)

    # invert Laplacian to obtain η
    (; pres) = progn.surface.timesteps[1]
    ∇⁻²!(pres, div, model.spectral_transform)

    # add perturbation
    η_perturb = SpeedyTransforms.spectral!(u, η_perturb_grid, model.spectral_transform)
    pres .+= η_perturb
    SpeedyTransforms.spectral_truncation!(pres)
    return nothing
end

end # module
