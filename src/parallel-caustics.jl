"""
    par_evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)

Does the same as [`evaluate_mass_homotopy`](@ref), but in a parallel way.

See also: [`calc_crit_curves`](@ref).
"""
function par_evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)
    mass_homotopy = create_mass_homotopy(masses, positions, E, Λ)
    init1, init2 = find_corrections(masses, positions, δs, E, Λ)
    start_roots = find_start_roots(mass_homotopy, positions, init1, init2, δs, find_root)
    hom(start) = homotopize(start, δs, 1., mass_homotopy, rate, mass_lim_func, find_root)
    roots = @showprogress 1 "Evaluating mass homotopy..." pmap(hom, start_roots)
    return roots
end


"""
    evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root)

Does the same as [`evaluate_angle_homotopy`](@ref), but in a parallel way.

See also: [`calc_crit_curves`](@ref).
"""
function par_evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root=simple_newton)
    angle_homotopy = create_angle_homotopy(masses, positions, E, Λ)
    crit_curves = Vector{Vector{Complex{Float64}}}(undef, 0)
    lim_func(t0, s0, rate) = 2π / nsteps
    hom(root) = homotopize_and_remember(root, 0., 2π, angle_homotopy, rate, lim_func, find_root)
    crit_curves = @showprogress 1 "Evaluating homotopy..." pmap(hom, roots)
    return crit_curves
end


"""
    par_calc_crit_curves(stars::Vector{Star}; E, Λ, δs=1e-6, rate=0.25, nsteps=200, find_root=simple_newton)

Does the same as [`calc_crit_curves`](@ref), but in a parallel way.
"""
function par_calc_crit_curves(stars::Vector{Star}; E, Λ, δs=1e-6, rate=0.1, nsteps=500, find_root=simple_newton)

    masses = [star.mass for star in stars]
    positions = [star.pos for star in stars]
    roots = par_evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)
    duplicate_warning_roots(roots)
    crit_curves = par_evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root)
    duplicate_warning_crit_curves(crit_curves)
    return crit_curves
end