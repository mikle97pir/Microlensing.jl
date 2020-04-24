mutable struct CellNode <: AbstractCell
       
    children::Vector{CellNode}
    nstars::Int
    
    mass::Float64 
    pos::Complex{Float64} # barycenter of the cell
    multipoles::Vector{Complex{Float64}}
    
    size::Float64
    center::Complex{Float64} # geometric center of the cell

    is_sink::Bool
    only_star::Star # makes sense only if this node is a leaf
    
    is_visited::Bool # technical variable for propagating multipoles

    @doc """
        CellNode()

    Creates a self-referential sink. It is a placeholder for non-existing children.
    """
    function CellNode()
       x = new()
       x.children = [x, x, x, x]
       x.is_sink = true
       return x
    end

    @doc """
        CellNode(sink::CellNode, size::Float64, order::Int)

    Creates a node without children. Array `children` consists of four `sink`s. `order` is the higher multipole degree taken into account.
    """
    function CellNode(sink::CellNode, size::Float64, order::Int)
       x = new()
       x.children = [sink, sink, sink, sink]
       x.nstars = 0
        
       x.mass = 0.
       x.pos = 0.
       x.multipoles = zeros(Complex{Float64}, order+1)
        
       x.size = size
       x.center = 0.
        
       x.is_sink = false
       x.only_star = Star(0, 0.)
        
       x.is_visited = false
       return x
    end
       
end


mutable struct CellTree
       
    sink::CellNode
    root::CellNode
    order::Int

    @doc """
        CellTree(size::Float64, order::Int)

    Creates a tree composed only of a root. All `root`'s children are set to `sink`.
    """
    function CellTree(size::Float64, order::Int)
       x = new()
       x.sink = CellNode()
       x.root = CellNode(x.sink, size, order)
       x.order = order  
       return x
    end
       
end


function Base.show(io::IO, cell::CellNode)
    if cell.is_sink
        print(io, "sink")
    else
        print(io, "cell")
    end  
end
    

function Base.show(io::IO, ::MIME"text/plain", cell::CellNode)
    if cell.is_sink
        print(io, "sink")
    else
        println(io, cell.children, "\n")
        println(io, "nstars = ", cell.nstars)
        println(io, "mass = ", cell.mass)
        println(io, "pos = ", cell.pos)
        println(io, cell.multipoles, "\n")
        println(io, "size = ", cell.size)
        println(io, "center = ", cell.center, "\n")
        println(io, "is_sink = ", cell.is_sink)
        println(io, "only_star = ", cell.only_star)
    end    
end


function AbstractTrees.children(cell::CellNode)
    if cell.is_sink
        return []
    else
        return cell.children
    end
end


"""
    next_index(cell::CellNode, star::Star)

Returns the index of `cell`'s child containing given `star`.
"""
function next_index(cell::CellNode, star::Star)
    
    Δx = star.pos.re - cell.center.re
    Δy = star.pos.im - cell.center.im
    
    if     (Δx < 0) & (Δy > 0)
        return 1
    elseif (Δx > 0) & (Δy > 0)
        return 2
    elseif (Δx < 0) & (Δy < 0)
        return 3
    elseif (Δx > 0) & (Δy < 0)
        return 4
    end               
end


"""
    child_center(cell::CellNode, i::Int)

Returns the center of the `i`-th child of the `cell`.
"""
function child_center(cell::CellNode, i::Int)
    mask =  (-1+im, 1+im, -1-im, 1-im)
    return cell.center + mask[i] * cell.size / 4
end


"""
    update_cell!(cell::CellNode, star::Star)

Updates `mass` and `pos` of the `cell` when the `star` is added.
"""
function update_cell!(cell::CellNode, star::Star)
    cell.nstars += 1
    M = cell.mass
    cell.mass += star.mass
    cell.pos *= M / cell.mass
    cell.pos += star.mass * star.pos / cell.mass 
end


"""
    add_cell!(T::CellTree, cell::CellNode, star::Star, i::Int)

Adds an `i`-th child containing `star` to the `cell`.

The square corresponding to the child *must* contain `star.pos`!
"""
function add_cell!(T::CellTree, cell::CellNode, star::Star, i::Int)
    child = CellNode(T.sink, cell.size / 2, T.order)
    update_cell!(child, star)
    child.only_star = star
    child.center = child_center(cell, i)
    cell.children[i] = child
end


"""
    add_star!(T::CellTree, star::Star)

Adds `star` to tree `T`.

It goes down the tree and modifies the cells until it meets a leaf. It then lengthens the branch until the stars are in different children.
"""
function add_star!(T::CellTree, star::Star)
    
    cell = T.root
    
    while cell.nstars > 1
        update_cell!(cell, star)
        i = next_index(cell, star)
        
        if !cell.children[i].is_sink
            cell = cell.children[i]
        else
            add_cell!(T, cell, star, i)
            break
        end
            
    end
    
    if cell.nstars == 1
        update_cell!(cell, star)
        old_star = cell.only_star
        while true
            i = next_index(cell, star)
            j = next_index(cell, old_star)
            
            if i != j
                add_cell!(T, cell, star, i)
                add_cell!(T, cell, old_star, j)
                break
            end
            
            if i == j
                add_cell!(T, cell, star, i)
                update_cell!(cell.children[i], old_star)
                cell = cell.children[i]
            end
        end        
    end
    
    if cell.nstars == 0
        update_cell!(cell, star)
        cell.only_star = star
    end
    
end


"""
    update_cell_multipoles!(cell::CellNode, child::CellNode, 
                            degrees::Vector{Complex{Float64}}, 
                            order::Int, binomials::Matrix{Int})

Updates `cell.multipoles` taking into account already known `child.multipoles`.

# Arguments
- `cell::CellNode`: the cell being updated.
- `child::CellNode`: the child being taken into account.
- `degrees::Vector{Complex{Float64}}`: vector of length `order+1`. It's a temporary array for keeping degrees of `A-B` from `0` to `order`.
- `order::Int`: the degree of the highest multipole taken into account.
- `binomials::Matrix{Int}`: matrix of signature `(order, order)`. `binomials[i, j]` should be equal to `binomial(j-1, i-1)`. It is normally computed by the [`calc_binomials`](@ref) function. 
"""
function update_cell_multipoles!(cell::CellNode, child::CellNode, 
                                 degrees::Vector{Complex{Float64}}, 
                                 order::Int, binomials::Matrix{Int})

    A = conj(child.pos)
    B = conj(cell.pos)

    degrees[1] = 1

    for i in 1:order
        degrees[i+1] = degrees[i]*(A-B)
    end

    cell.multipoles[1] += child.multipoles[1]

    for k in 1:order
        cell.multipoles[k+1] += -child.multipoles[1]*degrees[k+1]/k

        for j in 1:k
            P = binomials[j, k] * child.multipoles[j+1] * degrees[k-j+1]
            cell.multipoles[k+1] += P
        end

    end

end


"""
    calc_multipoles!(T::CellTree)

Calculates multipoles for the whole tree. All the stars should be already added. The function propagates multipoles from the leaves to the root.
"""
function calc_multipoles!(T::CellTree)
    
    S = Stack{CellNode}()
    push!(S, T.root)
    
    degrees = zeros(Complex{Float64}, T.order + 1)
    binomials = calc_binomials(T.order)
    
    while !isempty(S)
        cell = first(S)
        if cell.nstars > 1
            if !cell.is_visited
                for child in cell.children
                    if !child.is_sink
                        push!(S, child)
                    end
                end
                cell.is_visited = true
            else
                for child in cell.children
                    if !child.is_sink
                        update_cell_multipoles!(cell, child, 
                                                degrees, 
                                                T.order, binomials)
                    end
                end
                pop!(S)
            end
        end
            
        if cell.nstars == 1
            cell.multipoles[1] = cell.mass
            pop!(S)
        end
            
    end

end


"""
    build_tree(stars::Vector{Star}; width::Float64, order::Int=6)

Creates `CellTree` from a given array of stars. The root is centered at zero. 

`order` represents the highest multipole degree taken into account.
"""
function build_tree(stars::Vector{Star}; width::Float64, order::Int=6)
    T = CellTree(width, order)
    for star in stars
        add_star!(T, star)
    end
    calc_multipoles!(T)
    return T
end

