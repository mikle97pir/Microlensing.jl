"Every `AbstractCell` should have at least `center` and `size`. It should be square."
abstract type AbstractCell; end


struct Cell <: AbstractCell
    center::Complex{Float64}
    size::Float64
    leftup::Complex{Float64} # complex coordinate of the upper left corner
end


Cell(center, size) = Cell(center, size, center + (im-1)*size/2)


struct Grid
    center::Complex{Float64}
    size::Float64
    ngrid::Int # ngrid^2 is the number of cells in the grid
    step::Float64
    leftup::Complex{Float64}
end


function Grid(size, ngrid, center=0)
    step = size/ngrid
    leftup = center + (im-1)*size/2
    return Grid(center, size, ngrid, step, leftup)
end


Grid(cell::AbstractCell, ngrid) = Grid(cell.size, ngrid, cell.center)


"""
    get_index_leftup(grid::Grid, i, j)

Returns the complex coordinate of the upper left corner of the cell from the `i`-th row and `j`-th column.
"""
function get_index_leftup(grid::Grid, i, j)
    return grid.leftup + (j-1)*grid.step - (i-1)*grid.step*im
end


"""
    get_index_center(grid::Grid, i, j)

Returns the complex coordinate of the center of the cell from the `i`-th row and `j`-th column.
"""
function get_index_center(grid::Grid, i, j)
    return get_index_leftup(grid, i, j) + (1-im)*grid.step/2
end


"""
    Base.getindex(grid::Grid, i::Int, j::Int)

Returns the cell from the `i`-th row and `j`-th column.
"""
function Base.getindex(grid::Grid, i::Int, j::Int)
    return Cell(get_index_center(grid, i, j), grid.step)
end


"""
    matrix_rep(grid, n=grid.ngrid, kind=[:center or :leftup])

Returns a matrix composed of the complex coordinates of the cells from `grid`. They can be either centers or upper left corners depending on `kind` argument.

`n` argument tells where to stop: signature of resulting matrix is `(n, n)`. 
"""
function matrix_rep(grid::Grid, n=grid.ngrid; kind=:center)
    matrix = zeros(Complex{Float64}, (n, n))
    for i in 1:n
        for j in 1:n
            if kind == :center
                matrix[i, j] = get_index_center(grid, i, j)
            elseif kind == :leftup
                matrix[i, j] = get_index_leftup(grid, i, j)
            end
        end
    end
    return matrix
end


"See also: [`matrix_rep`](@ref)"
function matrix_rep!(matrix, grid::Grid, n=grid.ngrid; kind=:center)
    for i in 1:n
        for j in 1:n
            if kind == :center
                matrix[i, j] = get_index_center(grid, i, j)
            elseif kind == :leftup
                matrix[i, j] = get_index_leftup(grid, i, j)
            end
        end
    end
end


"""
    dist(c1::AbstractCell, c2::AbstractCell)

Calculates the distance between two cells as closed subsets of metric space.
"""
function dist(c1::AbstractCell, c2::AbstractCell)
    x1 = c1.center.re
    y1 = c1.center.im
    x2 = c2.center.re
    y2 = c2.center.im
    Δx = abs(x2 - x1)
    Δy = abs(y2 - y1)
    av_size = (c1.size + c2.size)/2
    return max(Δx - av_size, Δy - av_size, 0)
end