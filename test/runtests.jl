using Distributed
addprocs()
@everywhere using Microlensing
using Test
using Random

nstars = 50
rad = 6.
E = 1.
Λ = 1.

Random.seed!(0); stars = @test_logs generate_stars_ell(nstars, rad)

crit_curves = @test_logs calc_crit_curves(stars, E=E, Λ=Λ)
caustics = @test_logs calc_caustics(stars, crit_curves, E=E, Λ=Λ)

par_crit_curves = @test_logs par_calc_crit_curves(stars, E=E, Λ=Λ)
@test crit_curves == par_crit_curves

domain = @test_logs RectGrid(
    width=24.,
    nrows = 512,
    ncols = 512
)

image = @test_logs RectGrid(
    width=15.,
    nrows = 512,
    ncols = 512
)

tree = @test_logs build_tree(stars, width = 2*rad)

problem = @test_logs NumMLProblem(
    T = tree, 
    nstars = nstars, 
    nshare = 4, 
    nint = 4, 
    E = E, 
    Λ = Λ
)

mag = @test_logs calc_mag(problem, domain, image)
par_mag = @test_logs par_calc_mag(problem, domain, image)
@test mag == par_mag

domain1 = @test_logs RectGrid(12., (512, 256), -6.0+0.0im)
domain2 = @test_logs RectGrid(12., (512, 256), 6.0+0.0im)
mag1 = @test_logs par_calc_mag(problem, domain1, image)
mag2 = @test_logs par_calc_mag(problem, domain2, image)
@test mag == mag1 + mag2

image3 = @test_logs RectGrid(7.5, (512, 256), -3.75+0.0im)
image4 = @test_logs RectGrid(7.5, (512, 256), 3.75+0.0im)
mag3 = @test_logs par_calc_mag(problem, domain, image3)
mag4 = @test_logs par_calc_mag(problem, domain, image4)
@test mag == hcat(mag3, mag4)
