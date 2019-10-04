"""
    calc_mag(P::NumMLProblem, domain::Cell, image::Cell)

Just computes the magnification map on range `r`.
"""
function range_calc_mag(r::UnitRange{Int}, P::NumMLProblem, 
                           domain::Cell, image::Cell)

    ngrid, nshare, nint = P.ngrid, P.nshare, P.nint
    E, Λ = P.E, P.Λ

    mag = zeros(Int, (P.resol, P.resol))

    s_sig = (nshare + 1, nshare + 1)
    si_sig = (nshare*nint, nshare*nint)

    stack = Stack{CellNode}()
    near_stars = Vector{Star}(undef, P.nstars)
    far_sums = zeros(Complex{Float64}, s_sig)

    real_fs = zeros(Float64, s_sig)
    imag_fs = zeros(Float64, s_sig)

    int_far_sums = zeros(Complex{Float64}, si_sig)
    int_near_sums = zeros(Complex{Float64}, si_sig)
    int_grid_mat = zeros(Complex{Float64}, si_sig)

    lense = zeros(Complex{Float64}, si_sig)

    domain_grid = Grid(domain, ngrid)
    image_grid = Grid(image, P.resol)

    for i in r
        for j in 1:ngrid
            cell = domain_grid[i, j]

            nnstars = calc_far_sums!(far_sums, cell, P, near_stars, stack)

            int_grid = Grid(cell, nshare*nint)
            real_fs .= real.(far_sums)
            imag_fs .= imag.(far_sums)

            interpolate_far_sums!(int_far_sums, P, real_fs, imag_fs)

            matrix_rep!(int_grid_mat, int_grid, kind=:leftup)

            calc_near_sums!(int_near_sums, near_stars, 
                            nnstars, P, int_grid_mat)

            @. lense = E*Λ*real(int_grid_mat) + E*imag(int_grid_mat)*im
            @. lense = lense - int_far_sums - int_near_sums

            update_mag!(mag, lense, image_grid, P)

            fill!(far_sums, 0)
            fill!(int_near_sums, 0)
        end
    end
    return mag
end