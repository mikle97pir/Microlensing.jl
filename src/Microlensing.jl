module Microlensing

    import AbstractTrees
    using AbstractTrees: print_tree
    using DataStructures: Stack, top, pop!, push!
    using Parameters
    using Interpolations
    using Distributions
    using Random
    using LinearAlgebra: dot
    using ProgressMeter

    #### for distributions.jl
    using Distributions: @check_args, @distr_support
    import Distributions: shape, partype, pdf, logpdf, cdf
    import Base: maximum, minimum, convert, rand
    import Statistics: mean, quantile, var
    import StatsBase: params

    export Star, Power, Cell, NumMLProblem
    export generate_stars_rect, generate_stars_ell, build_tree, calc_mag,  range_calc_mag, find_crit_curves

    include("star-field/star-field.jl")
    include("grids-and-cells.jl")
    include("cell-trees.jl")
    include("magnification.jl")
    include("parallel-magnification.jl")
    include("caustics.jl")

end # module