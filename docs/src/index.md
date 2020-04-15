# Sources 

```@meta
CurrentModule = Microlensing
```

## Caustics

```@docs
    simple_newton(f, ∂f, init)
    find_corrections(masses, positions, δs, E, Λ)
    find_start_roots(mass_homotopy, positions, init1, init2, δs, find_root)
    homotopy_step(t0, s0, homotopy, Δs, find_root)
    auto_homotopy_step(t0, s0, homotopy, rate, lim_func, find_root)
    homotopize(t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)
    homotopize_and_remember!(t_list, s_list, t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)
    homotopize_and_remember(t_start, s_start, s_finish, homotopy, rate, lim_func, find_root)
    create_mass_homotopy(masses, positions, E, Λ)
    create_angle_homotopy(masses, positions, E, Λ)
    mass_lim_func(t0, s0, rate)
    evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)
    evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root)
    check_for_duplicates(array)
    duplicate_warning_crit_curves(crit_curves)
    duplicate_warning_roots(roots)
    calc_crit_curves(masses, positions, E, Λ, δs=1e-6, rate=0.25, nsteps=200, find_root=simple_newton)
    calc_caustics(masses, positions, E, Λ, crit_curves)
```

## Parallel caustics

```@docs
    par_evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)
    par_evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root=simple_newton)
    par_calc_crit_curves(masses, positions, E, Λ, δs=1e-6, rate=0.25, nsteps=200, find_root=simple_newton)
```

## Grids and cells

```@docs
    get_index_leftup(grid::Grid, i, j)
    get_index_center(grid::Grid, i, j)
    Base.getindex(grid::Grid, i::Int, j::Int)
    matrix_rep(grid::Grid, n=grid.ngrid; kind=:center)
    dist(c1::AbstractCell, c2::AbstractCell)
```

## Cell trees

```@docs
    CellNode()
    CellNode(sink::CellNode, size::Float64, order::Int)
    CellTree(size::Float64, order::Int)
    next_index(cell::CellNode, star::Star)
    child_center(cell::CellNode, i::Int)
    update_cell!(cell::CellNode, star::Star)
    add_cell!(T::CellTree, cell::CellNode, star::Star, i::Int)
    add_star!(T::CellTree, star::Star)
    update_cell_multipoles!(cell::CellNode, child::CellNode, degrees::Vector{Complex{Float64}}, order::Int, binomials::Matrix{Int})
    calc_binomials(order)
    calc_multipoles!(T::CellTree)
    build_tree(stars, size, order=6)
```

## Magnification

```@docs
    NumMLProblem
    calc_far_sums!(far_sums, cell::Cell, P::NumMLProblem, near_stars::Vector{Star}, stack::Stack{CellNode})
    interpolate_far_sums!(int_far_sums, P::NumMLProblem, real_fs, imag_fs)
    calc_near_sums!(int_near_sums, near_stars, nnstars, P::NumMLProblem, int_grid_mat)
    find_cell(pos, imsize, nimgrid)
    update_mag!(mag, lense, image_grid::Grid, P::NumMLProblem)
    normalize_mag(mag, P::NumMLProblem, domain::Cell, image::Cell)
    calc_mag(P::NumMLProblem, domain::Cell, image::Cell)
```

## Parallel magnification

```@docs
    range_calc_mag(r::UnitRange{Int}, P::NumMLProblem, domain::Cell, image::Cell)
    break_into_ranges(n, nranges)
    par_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)
    range_calc_mag!(mag, worker, r::UnitRange{Int}, P::NumMLProblem, 
                           domain::Cell, image::Cell)
    shared_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)
```

## Distributions

```@docs
    Power
```