module TemplateMatching

using StatsBase
using LinearAlgebra

export crosscorrelate

"""
    crosscorrelate(series, template, cc_eltype=Float64; normalize_template=false) -> Vector{cc_eltype}

    Return the normalized cross-correlation `cc` between `series` and `template`. If `series` is a vector of ``n`` elements 
    and template a vector of ``k`` elements, then `cc` is a vector of ``n - k + 1`` elements of type `cc_eltype`. 
    The element of `cc` at position ``l`` is Pearson correlation coefficient between the segment of `series` 
    starting at ``l`` and ``k`` elements long. See the [Wikipedia page](https://en.wikipedia.org/wiki/Cross-correlation#Normalization)
    on cross-correlation for more details.
"""
function crosscorrelate(series::AbstractVector{T1}, template::AbstractVector{T2}, 
    cc_eltype::Type{T3}=Float64; normalize_template=true) where {T1 <: Number, T2 <: Number, T3 <: AbstractFloat}
    
    if isempty(series) || isempty(template)
        error("Arguments must be non empty vectors")
    end

    if normalize_template
        y_mean, y_std = mean_and_std(template, corrected=false)
        y_norm = @. (template - y_mean) / y_std
    else
        y_norm = template
    end

    N = size(template, 1)
    cc = Vector{cc_eltype}(undef, size(series, 1) - N + 1)
    x_norm = similar(y_norm) 
    for n in eachindex(cc)
        x_view = view(series, n:n + N - 1)
        x_mean, x_std = mean_and_std(x_view, corrected=false)
        @. x_norm = (x_view - x_mean) / x_std
        cc[n] = dot(x_norm, y_norm) / N
    end
    cc
end

end
