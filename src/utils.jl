"""
    calc_binomials(order)

Computes binomials. Returns matrix of signature `(order, order)`, `binomials[i, j]` should be equal to `binomial(j-1, i-1)`. It is used in the [`calc_multipoles!`](@ref) function
."""
function calc_binomials(order)
    binomials = zeros(Int, order, order)
    for i in 0:(order-1)
       for j in 0:i
           binomials[j+1, i+1] = binomial(i, j)
       end
    end
    return binomials
end


"""
    simple_newton(f, ∂f, init)

Finds the root of the complex function `f` near the initial value `init`. `∂f` is the derivative of `f`. The function is a straightforward implementation of Newton method. It it used by default, for example, in [`calc_crit_curves`](@ref).
"""
function simple_newton(f, ∂f, init)
    T = init
    T_prev = init + 1.
    while T ≉ T_prev
        T_prev = T
        T = T_prev - f(T_prev)/∂f(T_prev)
        if isnan(T)
            error("The Newton method diverged. Try to decrease δs, decrease rate or increase nsteps. You may also try another root finder.")
        end
    end
    return T
end


"""
    check_for_duplicates(array)

Returns `true` if the `array` contains approximate(`≈`) duplicates. Else returnes `false`. It is used in the [`duplicate_warning_roots`](@ref) function.
"""
function check_for_duplicates(array)
    for i in 1:length(array)
        for j in (i+1):length(array)
            if array[i] ≈ array[j]
                return true
            end
        end
    end
    return false
end


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