"""
    range_calc_mag(r::UnitRange{Int}, P::NumMLProblem, domain::Cell, image::Cell)

Just computes the magnification map for the rows with numbers from the range `r`.
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



"""
    break_into_ranges(n, nranges)

Breaks the range `1:n` into a union of `nranges` ranges `1:x1`, `(x1+1):x2`, ..., `(xk+1):n`. Sizes of the resulting ranges cannot differ by more than 1. It is used by the [`par_calc_mag`](@ref) function to split the task between the workers.
"""
function break_into_ranges(n, nranges)
    x = div(n, nranges)
    y = x + 1
    a = n - nranges*x
    b = nranges - a
    aranges = [(1+(i-1)*y):(i*y) for i in 1:a]
    branges = [(a*y + 1 + (i-1)*x):(a*y + i*x) for i in 1:b]
    return vcat(aranges, branges)
end


"""
    par_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)

Does the same as [`calc_mag`](@ref), but in a parallel way. The progress bar does not work yet.
"""
function par_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)
    ngrid, nshare, nint = P.ngrid, P.nshare, P.nint
    E, Λ = P.E, P.Λ
    ranges = break_into_ranges(ngrid, nworkers())
    cmag(r) = range_calc_mag(r, P, domain, image)
    mag = sum(pmap(cmag, ranges))
    norm_mag = normalize_mag(mag, P, domain, image)
    return norm_mag
end


"""
    shared_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)

Does the same as [`par_calc_mag`](@ref), but uses a `SharedArray` for the result. It is a little bit slower, but uses memory efficiently.
"""
function shared_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)
    ngrid, nshare, nint = P.ngrid, P.nshare, P.nint
    E, Λ = P.E, P.Λ
    mag = SharedMatrix{Int}((P.resol, P.resol))
    ranges = break_into_ranges(ngrid, nworkers())
    progress_bar = Progress(ngrid, "Shooting rays...")
    channel = RemoteChannel(()->Channel{Bool}(ngrid), myid())
    @sync begin
        @async while true
            if take!(channel)
                next!(progress_bar)
                if progress_bar.counter == ngrid
                    break
                end
            end
        end
        for (i, worker) in enumerate(workers())
            @spawnat worker begin
                range_calc_mag!(mag, ranges[i], P, domain, image, channel)
                println("Finished")
            end
        end
    end
    norm_mag = normalize_mag(mag, P, domain, image)
    return norm_mag
end



"""
    range_calc_mag(r::UnitRange{Int}, P::NumMLProblem, domain::Cell, image::Cell)

Does the same as [`range_calc_mag`](@ref), but takes a `SharedArray` `mag` as an argument. `channel` argument is a `RemoteChannel` for updating the progress bar.

See [`shared_calc_mag`](@ref).
"""
function range_calc_mag!(mag, r::UnitRange{Int}, P::NumMLProblem, 
                           domain::Cell, image::Cell, channel)

    ngrid, nshare, nint = P.ngrid, P.nshare, P.nint
    E, Λ = P.E, P.Λ

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
        put!(channel, true)
    end
end