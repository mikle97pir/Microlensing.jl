"""
    find_corrections(masses, positions, δs, E, Λ)

Finds the approximate positions of the `2*nstars` roots of the jacobian, when all the masses are multiplied by a small number `δs`. The output arrays `init1` and `init2` should be used for the function [`find_start_roots`](@ref).

See also: [`calc_crit_curves`](@ref).
"""
function find_corrections(masses, positions, δs, E, Λ)
    inv_sums = [sum([1/(a - b) for b in positions if b != a]) for a in positions]
    U = δs*masses .* inv_sums/C(0, E, Λ)
    V = -δs*masses/(2*C(0, E, Λ))
    init1 = positions .- U .+ sqrt.(U.^2 .- 2*V) 
    init2 = positions .- U .- sqrt.(U.^2 .- 2*V)
    return init1, init2
end


"""
    find_start_roots(mass_homotopy, positions, init1, init2, δs, find_root)

Finds the exact positions of the `2*nstars` roots of the jacobian, when all the masses are multiplied by a small number `δs`. The `init1` and `init2` initial approximations should be computed by the function [`find_corrections`](@ref). 

`mass_homotopy` is a tuple of functions `(H, ∂tH, ∂sH, ∂ttH, ∂tsH)` containing the mass-changing homotopy `H(t, s)` and its derivatives. It should be generated by the function [`create_mass_homotopy`](@ref).

See also: [`calc_crit_curves`](@ref).
"""
function find_start_roots(mass_homotopy, positions, init1, init2, δs, find_root)
    H, ∂tH, ∂sH, ∂ttH, ∂tsH = mass_homotopy
    f(t) = H(t, δs)
    ∂f(t) = ∂tH(t, δs)
    nstars = length(positions)
    start_roots = zeros(Complex{Float64}, 2*nstars)
    @showprogress 1 "Searching for start roots..." for i in 1:nstars
        start_roots[2i-1] = find_root(f, ∂f, init1[i])
        start_roots[2i] = find_root(f, ∂f, init2[i])
    end
    return start_roots
end


"""
    homotopy_step(t0, s0, homotopy, Δs, find_root)

Finds `t` such that `F(t, s) = 0`, where `s = s0 + Δs` and `F(t0, s0) = 0`. Returns `t` and `s`.

`homotopy` is a tuple of functions `(F, ∂tF, ∂sF, ∂ttF, ∂tsF)` containing an arbitrary homotopy `F(t, s)` and its derivatives.

See also: [`calc_crit_curves`](@ref).
"""
function homotopy_step(t0, s0, homotopy, Δs, find_root)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = homotopy
    s = s0 + Δs
    Δt = -(∂sF(t0, s0)/∂tF(t0, s0)) * Δs
    f(t) = F(t, s)
    ∂f(t) = ∂tF(t, s)
    t = find_root(f, ∂f, t0+Δt)
    return t, s
end


"""
    auto_homotopy_step(t0, s0, homotopy, rate, lim_func, find_root)

Does the same as [`homotopy_step`](@ref), but the step `Δs` is determined automatically. 

`lim_func` is a function with arguments `t0, s0, rate` that limits the step from the above. For the mass-changing homotopy we use [`mass_lim_func`](@ref). For the angle-changing homotopy it is just a constant `2π / nsteps`.

See also: [`calc_crit_curves`](@ref).
"""
function auto_homotopy_step(t0, s0, homotopy, rate, lim_func, find_root)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = homotopy
    Δs = 2*rate*abs2(∂tF(t0, s0)) / abs(∂ttF(t0, s0)) / abs(∂sF(t0, s0))
     if Δs > lim_func(t0, s0, rate)
         Δs = lim_func(t0, s0, rate)
     end
    return homotopy_step(t0, s0, homotopy, Δs, find_root)
end


"""
    homotopize(t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)

Finds `t` such that `F(t, s_finish) = 0`, if it's known that `F(t_start, s_start) = 0`. Returns `t`. Does this by repeating [`auto_homotopy_step`](@ref) many times to get from `s_start` to `s_finish`. 

See also: [`calc_crit_curves`](@ref), [`homotopy_step`](@ref).
"""
function homotopize(t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = homotopy
    t0, s0 = t_start, s_start
    while s0 < s_finish
        t, s = auto_homotopy_step(t0, s0, homotopy, rate, lim_func, find_root)
        t0, s0 = t, s
    end
    Δs = s_finish - s0
    t, s = homotopy_step(t0, s0, homotopy, Δs, find_root)
    return t
end


"""
    homotopize_and_remember!(t_list, s_list, t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)

Does the same as [`homotopize`](@ref), but saves all the intermediate values of `t` and `s` to empty lists `t_list` and `s_list` respectively. 

See also: [`calc_crit_curves`](@ref), [`homotopy_step`](@ref).
"""
function homotopize_and_remember!(t_list, s_list, t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = homotopy
    t0, s0 = t_start, s_start
    while s0 < s_finish
        push!(t_list, t0)
        push!(s_list, s0)
        t, s = auto_homotopy_step(t0, s0, homotopy, rate, lim_func, find_root)
        t0, s0 = t, s
    end
    Δs = s_finish - s0
    t, s = homotopy_step(t0, s0, homotopy, Δs, find_root)
    push!(t_list, t)
    push!(s_list, s)
    return t
end


"""
    homotopize_and_remember(t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)

Does the same as [`homotopize_and_remember!`](@ref), but creates `t_list` and s_list` itself and returns `t_list`.

See also: [`calc_crit_curves`](@ref), [`homotopy_step`](@ref).
"""
function homotopize_and_remember(t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)
    t_list = Vector{Complex{Float64}}(undef, 0)
    s_list = Vector{Float64}(undef, 0)
    homotopize_and_remember!(t_list, s_list, t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)
    return t_list  
end


A(E, Λ) = E*(Λ + 1)/2
B(E, Λ) = E*(Λ - 1)/2
C(φ, E, Λ) = A(E, Λ)*exp(im*φ) - B(E, Λ)
∂C(φ, E, Λ) = im*A(E, Λ)*exp(im*φ)


"""
    create_mass_homotopy(masses, positions, E, Λ)

Returns the tuple `mass_homotopy` containing the mass-changing homotopy `H` and its derivatives `∂tH, ∂sH, ∂ttH, ∂tsH`.

See also: [`calc_crit_curves`](@ref).
"""
function create_mass_homotopy(masses, positions, E, Λ)
    H(t, s) = dot(s*masses, 1 ./ (t .- positions).^2) - C(0, E, Λ)
    ∂tH(t, s) = -2*dot(s*masses, 1 ./ (t .- positions).^3)
    ∂sH(t, s) = dot(masses, 1 ./ (t .- positions).^2)
    ∂ttH(t, s) = 6*dot(s*masses, 1 ./ (t .- positions).^4)
    ∂tsH(t, s) = -2*dot(masses, 1 ./ (t .- positions).^3)
    return H, ∂tH, ∂sH, ∂ttH, ∂tsH
end


"""
    create_angle_homotopy(masses, positions, E, Λ)

Returns the tuple `angle_homotopy` containing the angle-changing homotopy `G` and its derivatives `∂tG, ∂sG, ∂ttG, ∂tsG`.

See also: [`calc_crit_curves`](@ref).
"""
function create_angle_homotopy(masses, positions, E, Λ)
    G(t, s) = dot(masses, 1 ./ (t .- positions).^2) - C(s, E, Λ)
    ∂tG(t, s) = -2*dot(masses, 1 ./ (t .- positions).^3)
    ∂sG(t, s) = -∂C(s, E, Λ)
    ∂ttG(t, s) = 6*dot(masses, 1 ./ (t .- positions).^4)
    ∂tsG(t, s) = 0.
    return G, ∂tG, ∂sG, ∂ttG, ∂tsG
end

"""
    mass_lim_func(t0, s0, rate)

A limiting function for the step in the mass-changing homotopy.

See also: [`evaluate_mass_homotopy`](@ref), [`homotopize_and_remember!`](@ref), [`homotopize`](@ref), [`auto_homotopy_step`](@ref), [`calc_crit_curves`](@ref).
"""
mass_lim_func(t0, s0, rate) = (4/3)*rate*s0


"""
    evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)

Evaluates the mass-changing homotopy from δs to 1 and returns the array with the roots of the jacobian for the angle equal to zero.

See also: [`calc_crit_curves`](@ref).
"""
function evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)
    mass_homotopy = create_mass_homotopy(masses, positions, E, Λ)
    init1, init2 = find_corrections(masses, positions, δs, E, Λ)
    start_roots = find_start_roots(mass_homotopy, positions, init1, init2, δs, find_root)
    roots = similar(start_roots)
    @showprogress 1 "Evaluating mass homotopy..." for (i, start) in enumerate(start_roots)
        roots[i] = homotopize(start, δs, 1., mass_homotopy, rate, mass_lim_func, find_root)
    end
    return roots
end


"""
    evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root)

Evaluates the angle-changing homotopy from 0 to 2π and returns a `Vector{Vector{Complex{Float64}}}`, which consists of vectors of points lying on the critical curves.

See also: [`calc_crit_curves`](@ref).
"""
function evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root=simple_newton)
    angle_homotopy = create_angle_homotopy(masses, positions, E, Λ)
    crit_curves = Vector{Vector{Complex{Float64}}}(undef, 0)
    lim_func(t0, s0, rate) = 2π / nsteps
    @showprogress 1 "Evaluating homotopy..." for root in roots
        t_list = Vector{Complex{Float64}}(undef, 0)
        s_list = Vector{Float64}(undef, 0)
        homotopize_and_remember!(t_list, s_list, root, 0., 2π, angle_homotopy, rate, lim_func, find_root)
        push!(crit_curves, t_list)
    end
    return crit_curves
end


"""
    duplicate_warning_crit_curves(crit_curves)

Prints a warning if there are roots which glued to each other during the angle-changing homotopy.
"""
function duplicate_warning_crit_curves(crit_curves)
    last_points = [last(curve) for curve in crit_curves]
    if check_for_duplicates(last_points)
        @warn "Some roots stuck together, and the algorithm missed some critical curves. Try to decrease rate. You may also increase nsteps."
    end
    return nothing
end


"""
    duplicate_warning_roots(roots)

Prints a warning if there are roots which glued to each other during the mass-changing homotopy.
"""
function duplicate_warning_roots(roots)
    if check_for_duplicates(roots)
        @warn "Some roots stuck together. Try to decrease rate."
    end
    return nothing
end


"""
    calc_crit_curves(stars::Vector{Star}; E, Λ, δs=1e-6, rate=0.25, nsteps=200, find_root=simple_newton)

Normally returns a `Vector{Vector{Complex{Float64}}}`, which consists of vectors of points lying on the critical curves.

# Arguments
- `stars`: the array of stars.
- `E`, `Λ`: the lens parameters describing the external shear.
- `δs`: the starting parameter of the mass-changing homotopy; should be a small real number.
- `rate`: increasing this parameter effectively increases the adaptive step.
- `nsteps`: the minimal number of steps in the angle-changing homotopy; you should increase it if you want smoother curves.
- `find_root`: the root finder; its argumetns should be a complex function `f`, its derivative `∂f` and the initial approximation `init`.
"""
function calc_crit_curves(stars::Vector{Star}; E, Λ, δs=1e-6, rate=0.1, nsteps=500, find_root=simple_newton)

    masses = [star.mass for star in stars]
    positions = [star.pos for star in stars]
    roots = evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)
    duplicate_warning_roots(roots)
    crit_curves = evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root)
    duplicate_warning_crit_curves(crit_curves)
    return crit_curves
end


"""
    calc_caustics(stars::Vector{Star}, crit_curves; E, Λ)
Computes the caustics by given critical curves `crit_curves`. 

See also: [`calc_crit_curves`](@ref). 
"""
function calc_caustics(stars::Vector{Star}, crit_curves; E, Λ)
    caustics = Vector{Vector{Complex{Float64}}}(undef, 0)
    @showprogress 1 "Computing caustics..." for curve in crit_curves
        caustic = similar(curve)
        @showprogress 1 "Computing points on a caustic..." for (j, point) in enumerate(curve)
            res = A(E, Λ)*point + B(E, Λ)*conj(point)
            for star in stars
                res -= star.mass*(point-star.pos)/abs2(point-star.pos)
            end
            caustic[j] = res
        end
        push!(caustics, caustic)
    end
    return caustics
end

