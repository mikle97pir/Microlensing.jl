using Distributed
addprocs()
@everywhere using Microlensing
using Test
using Random

Random.seed!(0)
stars = @test_logs generate_stars_ell(50, 6)
masses = [star.mass for star in stars]
positions = [star.pos for star in stars]

crit_curves = @test_logs calc_crit_curves(masses, positions, 1., 1.)
caustics = @test_logs calc_caustics(masses, positions, 1., 1., crit_curves)

par_crit_curves = @test_logs par_calc_crit_curves(masses, positions, 1., 1.)
@test crit_curves == par_crit_curves

tree = @test_logs build_tree(stars, 24.)
problem = @test_logs NumMLProblem(
    T = tree, 
    nstars = 50, 
    nshare = 4, 
    nint = 4, 
    E = 1.0, 
    Λ = 1.0
)

domain = @test_logs RectGrid(24., (512, 512))
image = @test_logs RectGrid(15., (512, 512))

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