"""
    Power(α, a, b)

The *Power* distribution with shape `α` and support `[a, b]` has probability density function 

```math
    f(x, \\, \\alpha, \\, a, \\, b) = \\begin{cases}
        \\dfrac{1+\\alpha}{b^{1+\\alpha} - a^{1+\\alpha}} \\, x^{\\alpha}, \\quad \\alpha \\neq -1; \\\\
        \\dfrac{1}{\\log\\left(\\tfrac{b}{a}\\right) \\, x}, \\quad \\alpha = 1.
    \\end{cases}
```
"""
struct Power{T<:Real} <: ContinuousUnivariateDistribution
    α::T
    a::T
    b::T
    Power{T}(α::T, a::T, b::T) where {T} = new{T}(α, a, b)
end


function Power(α::T, a::T, b::T; check_args=true) where {T <: Real}
    check_args && @check_args(Power, a >= zero(a) && b > a)
    return Power{T}(α, a, b)
end


Power(α::Real, a::Real, b::Real) = Power(promote(α, a, b)...)
Power(α::Integer, a::Integer, b::Integer) = Power(float(α), float(a), float(b))


@distr_support Power d.a d.b


#### Conversions


convert(::Type{Power{T}}, α::Real, a::Real, b::Real) where {T<:Real} = Power(T(α), T(a), T(b))
convert(::Type{Power{T}}, d::Power{S}) where {T <: Real, S <: Real} = Power(T(d.α), T(d.a), T(d.b), check_args=false)


#### Parameters


shape(d::Power) = d.α
params(d::Power) = (d.α, d.a, d.b)
@inline partype(d::Power{T}) where {T<:Real} = T


#### Statistics


function mean(d::Power{T}) where T <: Real
    α, a, b = params(d)
    if α == -1
        return (b-a)/log(b/a)
    elseif α == -2
        return a*b*log(b/a)/(b-a)
    else
        return (a^(2+α)-b^(2+α))/(a^(1+α)-b^(1+α))*(1+α)/(2+α)
    end 
end


function var(d::Power{T}) where T <: Real
    α, a, b = params(d)
    if α == -1
        return ((a-b)*log(a/b)*(2*(a-b)+(a+b)*log(b/a)))/(2*log(b/a)^3)
    elseif α == -2
        return a*b*(1-(a*b*log(b/a)^2)/(a-b)^2)
    elseif α == -3
        return (2*a^2*b^2*(-2a+2b+(a+b)*log(a/b)))/((a-b)*(a+b)^2)
    else
        return ((1+α)*(a^(4+2α)+b^(4+2α)-a*b*(a*b)^α*(a^2*(2+α)^2+b^2*(2+α)^2- 2a*b*(1+α)*(3+α))))/((a^(1+α)-b^(1+α))^2*(2+α)^2*(3+α))
    end 
end  


#### Evaluation


function pdf(d::Power{T}, x::Real) where T<:Real
    α, a, b = params(d)
    if α == -1
        return 1/log(b/a) * 1/x
    else
        return (1+α)/(b^(1+α)-a^(1+α)) * x^α
    end
end


function logpdf(d::Power{T}, x::Real) where T<:Real
    α, a, b = params(d)
    if α == -1
        return -log(log(b/a)) - log(x)
    else
        return log((1+α)/(b^(1+α)-a^(1+α))) + α*log(x)
    end
end


function cdf(d::Power{T}, x::Real) where T<:Real
    α, a, b = params(d)
    if α == -1
        return log(x/a) / log(b/a)
    else
        return (a^(1+α) - x^(1+α)) / (a^(1+α) - b^(1+α))
    end
end


function quantile(d::Power{T}, q::Real) where T<:Real
    α, a, b = params(d)
    if α == -1
        return a*(b/a)^q
    else
        return (a^(1+α) - q*(a^(1+α) - b^(1+α))) ^ (1/(1+α))
    end
end


#### Sampling


function rand(rng::AbstractRNG, d::Power)
    α, a, b = params(d)
    q = (b-a)*rand(rng) + a
    return quantile(d, (q-a)/(b-a))
end


