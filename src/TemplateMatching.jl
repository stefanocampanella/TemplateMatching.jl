module TemplateMatching

using StatsBase
using LinearAlgebra

export crosscorrelate

"""
    crosscorrelate(series, template, cc_eltype=Float64; normalize_template=false)

Return the normalized cross-correlation between `series` and `template`. 

If `series` is a vector of ``n`` elements and template a vector of ``N`` elements, 
cross-correlation will be a vector of ``n - N + 1`` elements of type `cc_eltype`. 
The element of cross-correlation ``χ_k`` is Pearson correlation coefficient between 
``x^k_i`` denoting the elements of `series` ``x_{i + k}`` and 
the elements of `template` ``y_i``, with ``i = 1, 2, \\dots, N``

``
χ_k = \\frac{\\sum{\\left(x^k_i - μ_{x^k}\\right) \\left(y_i - μ_y\\right)}}{N σ_{x^k} σ_y} \\, ,
``

where ``μ`` and ``σ`` denotes the mean and variance respectively.

See the [Wikipedia page](https://en.wikipedia.org/wiki/Cross-correlation#Normalization) 
on cross-correlation for more details.

# Examples

```jldoctest
julia> crosscorrelate(sin.(0:0.25pi:2pi), [1, 1+√2, 1])
7-element Vector{Float64}:
  0.23258781949447394
  1.000000000000001
  0.23258781949447402
  7.401486830834377e-17
 -0.23258781949447394
 -1.0000000000000007
 -0.23258781949447394
```
"""
function crosscorrelate(series::AbstractVector{T1}, template::AbstractVector{T2}, 
    cc_eltype::Type{T3}=Float64; normalize_template=true) where {T1 <: Number, T2 <: Number, T3 <: AbstractFloat}
    
    if isempty(series) || isempty(template)
        error("Arguments must be non empty vectors")
    end

    if normalize_template
        template_mean, template_std = mean_and_std(template, corrected=false)
        y = @. (template - template_mean) / template_std
    else
        y = template
    end

    N = size(template, 1)
    cc = Vector{cc_eltype}(undef, size(series, 1) - N + 1)
    x = similar(y) 
    for n in eachindex(cc)
        series_view = view(series, n:n + N - 1)
        series_view_std = std(series_view, corrected=false)
        @. x = series_view / series_view_std
        cc[n] = dot(x, y) / N
    end
    cc
end

end
