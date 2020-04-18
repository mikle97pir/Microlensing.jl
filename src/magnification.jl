"""
# Fields
- `T::CellTree`: normally is created with the [`build_tree`](@ref) function.
- `nstars::Int`: number of stars in simulation.
- `ngrid::Int`: `ngrid^2` is the number of the first level cells. For them the "go down the tree and find the right nodes" part of the computation is excecuted completely.
- `nshare::Int`: `nshare^2` is the number of the second level subcells of a first level cell. For them the computation is complete, but they all share the same tree structure of the first level cell containing them.
- `nint::Int`: `nint^2` is the number of the third level subcells of a second level cell. For them the algorithm attentively takes care of the stars located near the corresponding first level cell, but the impact of all other stars is interpolated.
- `resol::Int`: `resol×resol` is the resolution of the resulting magnification map.
- `E::Float64`: first parameter of the shear.
- `Λ::Float64`: second parameter of the shear.
- `δ::Float64`: the precision parameter. Higher `δ` - higher precision!
"""
@with_kw struct NumMLProblem
    T::CellTree
    nstars::Int
    ngrid::Int
    nshare::Int = 8
    nint::Int = 8
    resol::Int
    E::Float64 = 1.
    Λ::Float64 = 1.
    δ::Float64 = 0.8
end


"""
    calc_far_sums!(far_sums, cell::Cell, P::NumMLProblem,
                        near_stars, stack)

Records to the `far_sums` matrix the impacts of far stars for all the second level subcells of the first level `cell`. Also saves all the stars located near the `cell` to the `near_stars` array.

`stack` is just a preallocated stack necessary for searching through the tree.
"""
function calc_far_sums!(far_sums, cell::Cell, P::NumMLProblem,
                        near_stars::Vector{Star}, stack::Stack{CellNode})

    nshare = P.nshare
    push!(stack, P.T.root)
    grid = Grid(cell, nshare)
    nnstars = 0

    while !isempty(stack)
        node = pop!(stack)

        if node.nstars == 1
            near_stars[nnstars + 1] = node.only_star
            nnstars += 1

            elseif dist(cell, node) / node.size > 1 / P.δ

                for i in 1:(nshare+1)
                    for j in 1:(nshare+1)
                        pos = get_index_leftup(grid, i, j)
                        vec = pos - node.pos
                        d2 = abs2(vec)
                        far_sums[i, j] += node.mass * vec / d2
                        power = vec/d2
                        for k in 2:(P.T.order+1)
                            power *= vec/d2
                            Δ = (k-1)*node.multipoles[k]*power
                            far_sums[i, j] -= Δ
                        end
                    end
                end

            else
                for child in node.children
                    if !child.is_sink
                        push!(stack, child)
                    end
                end
        end
    end

    return nnstars
end


"""
    interpolate_far_sums!(int_far_sums, P::NumMLProblem,
                               real_fs, imag_fs)

Calculates the interpolated impacts of the far stars for all the third level subcells of a first level cell and records them to the `int_far_sums` matrix. 

The `real_fs` and `imag_fs` matrices should contain correspondingly the real and imaginary parts of the far stars impact on the second level subcells. This impact is normally computed with the [`calc_far_sums!`](@ref) function.
"""
function interpolate_far_sums!(int_far_sums, P::NumMLProblem,
                               real_fs, imag_fs)

    nshare, nint = P.nshare, P.nint

    real_int = interpolate(real_fs, BSpline(Linear()))
    imag_int = interpolate(imag_fs, BSpline(Linear()))

    for i in 1:nint*nshare
        for j in 1:nint*nshare
            u = (i-1)/nint + 1
            v = (j-1)/nint + 1
            @inbounds Re = real_int(u, v)
            @inbounds Im = imag_int(u, v)
            @inbounds int_far_sums[i, j] = Re + Im*im
        end
    end

end


"""
    calc_near_sums!(int_near_sums, near_stars, 
                         nnstars, P::NumMLProblem, int_grid_mat)

Calculates the impact of all the stars close to a first level cell on its third level subcells. 

`nnstars` is the number of the near stars (it is less than `length(near_stars)`!) 

`int_grid_mat` is just a matrix of the upper left corners of the third level subcells.
"""
function calc_near_sums!(int_near_sums, near_stars, 
                         nnstars, P::NumMLProblem, int_grid_mat)

    nshare, nint = P.nshare, P.nint

    for i in 1:nint*nshare
        for j in 1:nint*nshare
            @inbounds for k in 1:nnstars
                point = int_grid_mat[i, j]
                star_pos = near_stars[k].pos
                r = point - star_pos 
                d2 = abs2(r)
                Δ = near_stars[k].mass * r/d2
                int_near_sums[i, j] += Δ
            end
        end
    end

end


"""
    find_cell(pos, imsize, nimgrid)

By the position `pos` of a point in the image plane finds the indicies of an image grid cell containing this point.
"""
function find_cell(pos, imsize, nimgrid)
    step = imsize/nimgrid
    x = pos.re
    y = pos.im
    i = ceil(Int, (imsize/2-y)/step)
    j = ceil(Int, (x+imsize/2)/step)
    return i, j
end


"""
    update_mag!(mag, lense, image_grid::Grid,
                     P::NumMLProblem)

Updates the magnification map `mag` taking into account the lense map `lense` computed for all the rays from a first level cell.
"""
function update_mag!(mag, lense, image_grid::Grid,
                     P::NumMLProblem)

    nshare, nint = P.nshare, P.nint

    for i in 1:nshare*nint
        for j in 1:nshare*nint
            pos = lense[i, j]
            imsize = image_grid.size
            nimgrid = image_grid.ngrid
            u, v = find_cell(pos, imsize, nimgrid)
            if (1 <= u <= nimgrid) & (1 <= v <= nimgrid)
                mag[u, v] += 1
            end
        end
    end

end


"""
    normalize_mag(mag, P::NumMLProblem, domain::Cell, image::Cell)
Normalizes the magnification map in such a way that magnification is equal to 1 on the infinity.
"""
function normalize_mag(mag, P::NumMLProblem, domain::Cell, image::Cell)
    nrays = P.ngrid^2*P.nshare^2*P.nint^2
    mult = nrays/(P.resol^2)*(image.size^2/domain.size^2)*(1/abs(P.Λ))
    return mag/mult
end


"""
    calc_mag(P::NumMLProblem, domain::Cell, image::Cell)

Just computes the magnification map. Output is normalized in such a way that the magnification is equal to 1 at the infinity.
"""
function calc_mag(P::NumMLProblem, domain::Cell, image::Cell)

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

    progress_bar = Progress(ngrid^2, "Shooting rays...")

    for i in 1:ngrid
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

            next!(progress_bar)
        end
    end
    
    return mag
end