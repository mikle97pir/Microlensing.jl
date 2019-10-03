include("distributions.jl")

struct Star
    mass::Float64
    pos::Complex{Float64}
end


function generate_stars_rect(nstars, a, b=a; dmass=DiscreteUniform(1, 1)) 
    x = a*rand(nstars) .- a/2
    y = b*rand(nstars) .- b/2
    m = rand(dmass, nstars)
    return Star.(m, x .+ y*im)
end


function generate_stars_ell(nstars, a, b=a; φ=0, 
                            dpos=Power(1, 0, 1), dmass=DiscreteUniform(1, 1)) 
    r = a*rand(dpos, nstars)
    θ = 2π*rand(nstars)
    x = r .* cos.(θ)
    y = r .* sin.(θ)
    y .= y .* b/a
    z = x .+ y*im
    z .= z .* exp(im*φ)
    m = rand(dmass, nstars)
    return Star.(m, z)
end