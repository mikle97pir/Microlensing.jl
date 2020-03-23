function find_root(f, ∂f, init, prec=1e-10)
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


function find_corrections(masses, coords, Δs, E, Λ)
    inv_sums = [sum([1/(a - b) for b in coords if b != a]) for a in coords]
    U = Δs*masses .* inv_sums/C(0, E, Λ)
    V = -Δs*masses/(2*C(0, E, Λ))
    init1 = coords .- U .+ sqrt.(U.^2 .- 2*V) 
    init2 = coords .- U .- sqrt.(U.^2 .- 2*V)
    return init1, init2
end


function find_start_roots(mass_FD, coords, init1, init2, Δs)
    H, ∂tH, ∂sH, ∂ttH, ∂tsH = mass_FD
    f(t) = H(t, Δs)
    ∂f(t) = ∂tH(t, Δs)
    nstars = length(coords)
    start_roots = zeros(Complex{Float64}, 2*nstars)
    for i in 1:nstars
        start_roots[2i-1] = find_root(f, ∂f, init1[i])
        start_roots[2i] = find_root(f, ∂f, init2[i])
    end
    return start_roots
end


function homotopy_step(t0, s0, FD, Δs)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = FD
    s = s0 + Δs
    Δt = -(∂sF(t0, s0)/∂tF(t0, s0)) * Δs
    f(t) = F(t, s)
    ∂f(t) = ∂tF(t, s)
    t = find_root(f, ∂f, t0+Δt)
    return t, s
end


function auto_homotopy_step(t0, s0, FD, rate, lim_func)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = FD
    Δs = rate*abs2(∂tF(t0, s0)) / abs(∂ttF(t0, s0)) / abs(∂sF(t0, s0))
     if Δs > lim_func(t0, s0, rate)
         Δs = lim_func(t0, s0, rate)
     end
    return homotopy_step(t0, s0, FD, Δs)
end


function homotopize(t_start, s_start, s_finish, FD, rate, lim_func)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = FD
    t0, s0 = t_start, s_start
    while s0 < s_finish
        t, s = auto_homotopy_step(t0, s0, FD, rate, lim_func)
        t0, s0 = t, s
    end
    Δs = s_finish - s0
    t, s = homotopy_step(t0, s0, FD, Δs)
    return t
end


function homotopize_and_remember!(t_list, s_list, t_start, s_start, s_finish, FD, rate, lim_func)
    F, ∂tF, ∂sF, ∂ttF, ∂tsF = FD
    t0, s0 = t_start, s_start
    while s0 < s_finish
        push!(t_list, t0)
        push!(s_list, s0)
        t, s = auto_homotopy_step(t0, s0, FD, rate, lim_func)
        t0, s0 = t, s
    end
    Δs = s_finish - s0
    t, s = homotopy_step(t0, s0, FD, Δs)
    push!(t_list, t)
    push!(s_list, s)
    return t
end


A(E, Λ) = E*(Λ + 1)/2
B(E, Λ) = E*(Λ - 1)/2
C(φ, E, Λ) = A(E, Λ)*exp(im*φ) - B(E, Λ)
∂C(φ, E, Λ) = im*A(E, Λ)*exp(im*φ)


function create_mass_homotopy(masses, coords, E, Λ)
    H(t, s) = dot(s*masses, 1 ./ (t .- coords).^2) - C(0, E, Λ)
    ∂tH(t, s) = -2*dot(s*masses, 1 ./ (t .- coords).^3)
    ∂sH(t, s) = dot(masses, 1 ./ (t .- coords).^2)
    ∂ttH(t, s) = 6*dot(s*masses, 1 ./ (t .- coords).^4)
    ∂tsH(t, s) = -2*dot(masses, 1 ./ (t .- coords).^3)
    return H, ∂tH, ∂sH, ∂ttH, ∂tsH
end


function create_homotopy(masses, coords, E, Λ)
    G(t, s) = dot(masses, 1 ./ (t .- coords).^2) - C(s, E, Λ)
    ∂tG(t, s) = -2*dot(masses, 1 ./ (t .- coords).^3)
    ∂sG(t, s) = -∂C(s, E, Λ)
    ∂ttG(t, s) = 6*dot(masses, 1 ./ (t .- coords).^4)
    ∂tsG(t, s) = 0.
    return G, ∂tG, ∂sG, ∂ttG, ∂tsG
end

mass_lim_func(t0, s0, rate) = (2/3)*rate*s0

function evaluate_mass_homotopy(masses, coords, E, Λ, Δs=1e-6, rate=0.5)
    mass_FD = create_mass_homotopy(masses, coords, E, Λ)
    init1, init2 = find_corrections(masses, coords, Δs, E, Λ)
    start_roots = find_start_roots(mass_FD, coords, init1, init2, Δs)
    roots = similar(start_roots)
    @showprogress 1 "Computing..." for (i, start) in enumerate(start_roots)
        roots[i] = homotopize(start, Δs, 1., mass_FD, rate, mass_lim_func)
    end
    return roots
end


function evaluate_homotopy(roots, masses, coords, E, Λ, nsteps=100, rate=0.5)
    FD = create_homotopy(masses, coords, E, Λ)
    crit_curves = Vector{Vector{Complex{Float64}}}(undef, 0)
    lim_func(t0, s0, rate) = 2π / nsteps
    @showprogress 1 "Computing..." for root in roots
        t_list = Vector{Complex{Float64}}(undef, 0)
        s_list = Vector{Float64}(undef, 0)
        homotopize_and_remember!(t_list, s_list, root, 0., 2π, FD, 0.1, lim_func)
        push!(crit_curves, t_list)
    end
    return crit_curves
end


function check_mass_homotopy(roots)
    for i in 1:length(roots)
        for j in (i+1):length(roots)
            if roots[i]≈roots[j]
                println("bad pair of roots: ")
                @show (i, j)
            end
        end
    end
end