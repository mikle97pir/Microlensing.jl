using Microlensing
using Test
using Random

Random.seed!(0)
stars = @test_logs generate_stars_ell(50, 6)
masses = [star.mass for star in stars]
positions = [star.pos for star in stars]

crit_curves = @test_logs calc_crit_curves(masses, positions, 1., 1.)
caustics = @test_logs calc_caustics(masses, positions, 1., 1., crit_curves)

tree = @test_logs build_tree(stars, 24.)
problem = @test_logs NumMLProblem(tree, 50, 512, 4, 4, 512, 1., 1., 0.8)
domain = @test_logs Cell(0., 24.)
image = @test_logs Cell(0., 15.)

mag = @test_logs calc_mag(problem, domain, image)