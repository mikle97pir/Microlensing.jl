"""
    break_into_ranges(n::Int, nranges::Int)

Breaks the range `1:n` into a union of `nranges` ranges `1:x1`, `(x1+1):x2`, ..., `(xk+1):n`. Sizes of the resulting ranges cannot differ by more than 1. It is used by the [`par_calc_mag`](@ref) function to split the task between the workers.
"""
function break_into_ranges(n::Int, nranges::Int)
    x = div(n, nranges)
    y = x + 1
    a = n - nranges*x
    b = nranges - a
    aranges = [(1+(i-1)*y):(i*y) for i in 1:a]
    branges = [(a*y + 1 + (i-1)*x):(a*y + i*x) for i in 1:b]
    return vcat(aranges, branges)
end


"""
    range_calc_mag!(mag, r::UnitRange{Int}, P::NumMLProblem, domain::Cell, image::Cell, channel::RemoteChannel{Channel{Bool}})

Computes the magnification map for a part of the set of rays (with row numbers from the range `r`) and adds it to `mag`. It is used by [`par_calc_mag`](@ref) to do a part of the job on one of the workers. The `channel` is used to transmit information about the calculation progress. 
"""
function range_calc_mag!(mag, r::UnitRange{Int}, P::NumMLProblem, domain::Cell, image::Cell, channel::RemoteChannel{Channel{Bool}})

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


"""
    par_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)

Does the same as [`calc_mag`](@ref), but in a parallel way.
"""
function par_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)
    ngrid, nshare, nint = P.ngrid, P.nshare, P.nint
    E, Λ = P.E, P.Λ
    mag = SharedMatrix{Int}((P.resol, P.resol))
    ranges = break_into_ranges(ngrid, nworkers())
    progress_bar = Progress(ngrid, "Shooting rays...")
    channel = RemoteChannel(()->Channel{Bool}(ngrid), myid())
    @sync begin
        # this async block updates the progress_bar
        @async while true
            if take!(channel)
                next!(progress_bar)
                if progress_bar.counter == ngrid
                    break
                end
            end
        end
        # here starts the real computation
        for (i, worker) in enumerate(workers())
            @spawnat worker begin
                range_calc_mag!(mag, ranges[i], P, domain, image, channel)
            end
        end
    end
    return sdata(mag)
end



"""
    distr_calc_mag(P::NumMLProblem, domain::Cell, image::Cell)

Does the same as [`calc_mag`](@ref), but in a parallel way.
"""
function distr_calc_mag(P::NumMLProblem, domain::Cell, image::Cell, tmp_path="./")
    ngrid, nshare, nint = P.ngrid, P.nshare, P.nint
    E, Λ = P.E, P.Λ
    mag = dzeros(
        Int, (P.resol, P.resol, nworkers()), 
        workers(), [1, 1, nworkers()]
    )
    ranges = break_into_ranges(ngrid, nworkers())
    progress_bar = Progress(ngrid, "Shooting rays...")
    channel = RemoteChannel(()->Channel{Bool}(ngrid), myid())

    @sync begin
        # this async block updates the progress_bar
        @async while true
            if take!(channel)
                next!(progress_bar)
                if progress_bar.counter == ngrid
                    break
                end
            end
        end
        # here starts the real computation
        for (i, worker) in enumerate(workers())
            @spawnat worker begin
                range_calc_mag!(mag, ranges[i], P, domain, image, channel)
            end
        end
    end

    @info "Computation finished"
    sleep(1)

    # saving the magnification maps to temporary files
    temp_save_mags(mag, tmp_path)
    close(mag)
    @info "Temporary files created"
    sleep(1)

    #cleaning up
    @everywhere GC.gc()

    @info "Memory is free"

    sleep(1)

    # loading and adding up all the maps from workers
    return sum_mags(P.resol, tmp_path)
end


function temp_save_mags(mag::DArray, tmp_path="./")
    @sync for worker in workers()
        @spawnat worker begin
            path = tmp_path*"tmp"*string(worker)*".jld"
            jldopen(path, "w") do file
                write(file, "mag", view(mag.localpart, :, :, 1))
            end
        end
    end
end

function sum_mags(resol, tmp_path="./")
    mag = zeros(Int, (resol, resol))
    @showprogress 1 "Fetching from workers..." for worker in workers()
        path = tmp_path*"tmp"*string(worker)*".jld"
        jldopen(path, "r") do file
            mag .+= read(file, "mag")
        end
        rm(path)
    end
    return mag
end