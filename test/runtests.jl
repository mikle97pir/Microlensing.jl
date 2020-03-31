using Microlensing
using Test
using Random

Random.seed!(0)
stars = @test_logs generate_stars_ell(200, 15)
masses = [star.mass for star in stars]
positions = [star.pos for star in stars]

crit_curves = @test_logs calc_crit_curves(masses, positions, 1., 1.)
caustics = @test_logs calc_caustics(masses, positions, 1., 1., crit_curves)

tree = @test_logs build_tree(stars, 60.)
problem = @test_logs NumMLProblem(tree, 200, 1024, 4, 4, 1024, 1., 1., 0.8)
domain = @test_logs Cell(0., 60)
image = @test_logs Cell(0., 15.)

mag = @test_logs calc_mag(problem, domain, image)