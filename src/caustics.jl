function find_root(f, ∂f, init, prec=100*eps())
    ΔT = prec + 1
    T = init
    n_iter = 0
    while ΔT > prec
        T_next = T - f(T)/∂f(T)
        ΔT = abs(T_next - T)
        T = T_next
        n_iter += 1
    end
    if isnan(T)
        throw(ErrorException)
    end
    return T
end


function grid_homotopize!(result, G, ∂tG, ∂sG, start, finish, ngrid)
    step = (finish - start)/ngrid
    for i in 2:(ngrid+1)
        s = start + (i-1)*step
        f(t) = G(t, s)
        ∂f(t) = ∂tG(t, s)
        for (j, root) in enumerate(result[i-1, :])
            init = root - ∂sG(root, s-step)/∂tG(root,s-step) * step
            result[i, j] = find_root(f, ∂f, init)
        end
    end
end


function adaptive_homotopize!(roots, next_roots, H, ∂tH, ∂sH, 
                              start, finish, init_step)
    k = 0
    step = init_step
    s = start
    while true
        if k == 10
            step *= 5
            k = 0
        end
        try
            if s >= finish
                step = finish - s
            end
            s_next = s + step
            f(t) = H(t, s_next)
            ∂f(t) = ∂tH(t, s_next)
            for (i, root) in enumerate(roots)
                init = root - ∂sH(root, s)/∂tH(root, s) * step
                next_roots[i] = find_root(f, ∂f, init)
            end
            k +=1
            roots .= next_roots
            s = s_next
            if s ≈ finish
                return s
            end
        catch ErrorException
            step /= 10
            k = 0
        end
    end
end


function find_corrections(masses, coords, s, E, Λ)
    inv_sums = [sum([1/(a - b) for b in coords if b != a]) for a in coords]
    U = s*masses .* inv_sums/C(0, E, Λ)
    V = -s*masses/(2*C(0, E, Λ))
    corr1 = -U .+ sqrt.(U.^2 .- 2*V) 
    corr2 = -U .- sqrt.(U.^2 .- 2*V) 
    return corr1, corr2
end


function find_start_roots!(roots, H, ∂tH, start, nstars, coords, corr1, corr2)
    f(t) = H(t, start)
    ∂f(t) = ∂tH(t, start)
    for i in 1:nstars
        init1 = coords[i]+corr1[i]
        init2 = coords[i]+corr2[i]
        roots[2i-1] = find_root(f, ∂f, init1)
        roots[2i] = find_root(f, ∂f, init2)
    end
end


A(E, Λ) = E*(Λ + 1)/2
B(E, Λ) = E*(Λ - 1)/2
C(φ, E, Λ) = A(E, Λ)*exp(im*φ) - B(E, Λ)
∂C(φ, E, Λ) = im*A(E, Λ)*exp(im*φ)


function create_mass_homotopy(masses, coords, E, Λ)
    H(t, s) = dot(s*masses, 1 ./ (t .- coords).^2) - C(0, E, Λ)
    ∂tH(t, s) = -2*dot(s*masses, 1 ./ (t .- coords).^3)
    ∂sH(t, s) = dot(masses, 1 ./ (t .- coords).^2)
    return H, ∂tH, ∂sH
end


function create_homotopy(masses, coords, E, Λ)
    G(t, s) = dot(masses, 1 ./ (t .- coords).^2) - C(s, E, Λ)
    ∂tG(t, s) = -2*dot(masses, 1 ./ (t .- coords).^3)
    ∂sG(t, s) = ∂C(s, E, Λ)
    return G, ∂tG, ∂sG
end


function find_mass_start(nstars, masses, coords, E, Λ)
    min_dist2 = abs2(coords[2]-coords[1])
    max_mass = masses[1]
    for i in 1:nstars
        if masses[i] > max_mass
            max_mass = masses[i]
        end
        for j in (i+1):nstars
            dist2 = abs2(coords[i]-coords[j])
            if dist2 < min_dist2
                min_dist2 = dist2
            end
        end
    end
    mass_bound = min_dist2/nstars*abs(C(0, E, Λ))
    return 0.1*mass_bound/max_mass
end


function calc_crit_curves(nstars, stars, E, Λ, nφ)
    masses = [star.mass for star in stars]
    coords = [star.pos for star in stars]
    roots = zeros(Complex{Float64}, 2*nstars)
    next_roots = zeros(Complex{Float64}, 2*nstars)
    crit_curves = zeros(Complex{Float64}, (nφ+1, 2*nstars))

    H, ∂tH, ∂sH = create_mass_homotopy(masses, coords, E, Λ)
    G, ∂tG, ∂sG = create_homotopy(masses, coords, E, Λ)

    start = find_mass_start(nstars, masses, coords, E, Λ)
    corr1, corr2 = find_corrections(masses, coords, start, E, Λ)
    find_start_roots!(roots, H, ∂tH, start, nstars, coords, corr1, corr2)
    adaptive_homotopize!(roots, next_roots, H, ∂tH, ∂sH, start, 1., start)

    crit_curves[1, :] = roots 
    grid_homotopize!(crit_curves, G, ∂tG, ∂sG, 0, 2π, nφ)
    return crit_curves
end


function calc_caustics(nstars, stars, nφ, E, Λ, crit_curves)
    caustics = similar(crit_curves)
    for i in 1:2*nstars
        for j in 1:(nφ+1)
            point = crit_curves[j, i]
            res = A(E, Λ)*point + B(E, Λ)*conj(point)
            for star in stars
                res -= star.mass*(point-star.pos)/abs2(point-star.pos)
            end
            caustics[j, i] = res
        end
    end
    return caustics
end