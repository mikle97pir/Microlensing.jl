# Sources 

```@meta
CurrentModule = Microlensing
```
## Utils

```@docs
    calc_binomials(order)
    simple_newton(f, ∂f, init)
    check_for_duplicates(array)
    break_into_ranges(n::Int, nranges::Int)
```

## Caustics

```@docs
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
    duplicate_warning_crit_curves(crit_curves)
    duplicate_warning_roots(roots)
    calc_crit_curves(stars::Vector{Star}; E, Λ, δs=1e-6, rate=0.1, nsteps=500, find_root=simple_newton)
    calc_caustics(stars::Vector{Star}, crit_curves; E, Λ)
```

## Parallel caustics

```@docs
    par_evaluate_mass_homotopy(masses, positions, E, Λ, δs, rate, find_root)
    par_evaluate_angle_homotopy(roots, masses, positions, E, Λ, rate, nsteps, find_root=simple_newton)
    par_calc_crit_curves(stars::Vector{Star}; E, Λ, δs=1e-6, rate=0.1, nsteps=500, find_root=simple_newton)
```

## Grids and cells

```@docs
    get_index_leftup(grid::AbstractGrid, i, j)
    get_index_center(grid::AbstractGrid, i, j)
    Base.getindex(grid::AbstractGrid, i::Int, j::Int)
    matrix_rep(grid::RectGrid, n=grid.ngrid; kind=:center)
    matrix_rep!(matrix, grid::RectGrid, n=grid.ngrid; kind=:center)
    matrix_rep(grid::SquareGrid, n=grid.ngrid; kind=:center)
    matrix_rep!(matrix, grid::SquareGrid, n=grid.ngrid; kind=:center)
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
    calc_multipoles!(T::CellTree)
    build_tree(stars::Vector{Star}; width::Float64, order::Int=6)
```

## Magnification

```@docs
    NumMLProblem
    calc_far_sums!(far_sums, cell::Cell, P::NumMLProblem, near_stars::Vector{Star}, stack::Stack{CellNode})
    interpolate_far_sums!(int_far_sums, P::NumMLProblem, real_fs, imag_fs)
    calc_near_sums!(int_near_sums, near_stars, nnstars, P::NumMLProblem, int_grid_mat)
    find_cell(pos, image::RectGrid)
    update_mag!(mag, lense, image::RectGrid, P::NumMLProblem)
    update_mag!(mag::DArray, lense, image::RectGrid, P::NumMLProblem)
    calc_mag(P::NumMLProblem, domain::RectGrid, image::RectGrid)
```

## Parallel magnification

```@docs
    range_calc_mag!(mag, r::UnitRange{Int}, P::NumMLProblem, domain::RectGrid, image::RectGrid, channel::RemoteChannel{Channel{Bool}})
    shared_calc_mag(P::NumMLProblem, domain::RectGrid, image::RectGrid)
    par_calc_mag(P::NumMLProblem, domain::RectGrid, image::RectGrid, tmp_path="./")
    temp_save_mags(mag::DArray, tmp_path="./")
    sum_mags(image::RectGrid, tmp_path="./")
```

## Distributions

```@docs
    Power
```