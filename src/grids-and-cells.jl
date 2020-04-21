"Every `AbstractCell` should have at least `center` and `size`. It should be square."
abstract type AbstractCell; end


struct Cell <: AbstractCell
    center::Complex{Float64}
    size::Float64
    leftup::Complex{Float64} # complex coordinate of the upper left corner
end


Cell(center, size) = Cell(center, size, center + (im-1)*size/2)


abstract type AbstractGrid; end

struct RectGrid <: AbstractGrid
    center::Complex{Float64}
    size::Tuple{Float64, Float64}
    ngrid::Tuple{Int, Int} 
    # ngrid[1]*ngrid[2] is the number of cells in the grid
    step::Float64
    leftup::Complex{Float64}
end


function RectGrid(width::Float64, ngrid::Tuple{Int, Int}, 
                                  center::Complex{Float64}=0.0+0.0im)
    step = width / ngrid[2]
    height = step * ngrid[1]
    size = (width, height)
    leftup = center + (im*height - width)/2
    return RectGrid(center, size, ngrid, step, leftup)
end


struct SquareGrid <: AbstractGrid
    center::Complex{Float64}
    size::Float64
    ngrid::Int # ngrid^2 is the number of cells in the grid
    step::Float64
    leftup::Complex{Float64}
end


function SquareGrid(size::Float64, ngrid::Int, 
                    center::Complex{Float64}=0.0+0.0im)
    step = size/ngrid
    leftup = center + (im-1)*size/2
    return SquareGrid(center, size, ngrid, step, leftup)
end


SquareGrid(cell::Cell, ngrid::Int) = SquareGrid(cell.size, ngrid, cell.center)


"""
    get_index_leftup(grid::SquareGrid, i, j)

Returns the complex coordinate of the upper left corner of the cell from the `i`-th row and `j`-th column.
"""
function get_index_leftup(grid::AbstractGrid, i, j)
    return grid.leftup + (j-1)*grid.step - (i-1)*grid.step*im
end


"""
    get_index_center(grid::AbstractGrid, i, j)

Returns the complex coordinate of the center of the cell from the `i`-th row and `j`-th column.
"""
function get_index_center(grid::AbstractGrid, i, j)
    return get_index_leftup(grid, i, j) + (1-im)*grid.step/2
end


"""
    Base.getindex(grid::AbstractGrid, i::Int, j::Int)

Returns the cell from the `i`-th row and `j`-th column.
"""
function Base.getindex(grid::AbstractGrid, i::Int, j::Int)
    return Cell(get_index_center(grid, i, j), grid.step)
end


"""
    matrix_rep(grid::RectGrid, n=grid.ngrid, kind=[:center or :leftup])

See also: [`matrix_rep!`](@ref)
"""
function matrix_rep(grid::RectGrid, n=grid.ngrid; kind=:center)
    matrix = zeros(Complex{Float64}, n)
    matrix_rep!(matrix, grid, n, kind=kind)
    return matrix
end


"""
    matrix_rep!(matrix, grid::RectGrid, n=grid.ngrid; kind=[:center or :leftup])

Fills the `matrix` with the complex coordinates of the cells from `grid`. They can be either centers or upper left corners depending on `kind` argument.

`n` argument tells where to stop: signature of resulting matrix is equal to `n`. 
"""
function matrix_rep!(matrix, grid::RectGrid, n=grid.ngrid; kind=:center)
    for i in 1:n[1]
        for j in 1:n[2]
            if kind == :center
                matrix[i, j] = get_index_center(grid, i, j)
            elseif kind == :leftup
                matrix[i, j] = get_index_leftup(grid, i, j)
            end
        end
    end
end


"""
    matrix_rep(grid::SquareGrid, n=grid.ngrid, kind=[:center or :leftup])

See also: [`matrix_rep!`](@ref)
"""
function matrix_rep(grid::SquareGrid, n=grid.ngrid; kind=:center)
    matrix = zeros(Complex{Float64}, (n, n))
    matrix_rep!(matrix, grid, n, kind=kind)
    return matrix
end


"""
    matrix_rep!(matrix, grid::SquareGrid, n=grid.ngrid; kind=[:center or :leftup])

Fills the `matrix` with the complex coordinates of the cells from `grid`. They can be either centers or upper left corners depending on `kind` argument.

`n` argument tells where to stop: signature of resulting matrix is `(n, n)`. 
"""
function matrix_rep!(matrix, grid::SquareGrid, n=grid.ngrid; kind=:center)
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