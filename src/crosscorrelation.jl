"""
    crosscorrelate(series, template, [element_type=Float64], normalize_template=true)

Return the normalized cross-correlation between `series` and `template`. 

If `series` is a vector of ``n`` elements and template a vector of ``N`` elements, 
cross-correlation will be a vector of ``n - N + 1`` elements of type `element_type`. 
The element of cross-correlation ``χ_k`` is Pearson correlation coefficient between 
``x^{(k)}_i`` denoting the elements of `series` ``x_{i + k}`` and 
the elements of `template` ``y_i``, with ``i = 1, 2, \\dots, N``

``
χ_k = \\frac{\\sum{\\left(x^{(k)}_i - μ_{x^{(k)}}\\right) 
\\left(y_i - μ_y\\right)}}{N σ_{x^{(k)}} σ_y} \\, ,
``

where ``μ`` and ``σ`` denotes the mean and variance respectively. 

If `normalize_template` is set to false then `template` is assumed to have mean and std
respectively equal to 0.0 and 1.0.

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
function crosscorrelate(series::AbstractVector{T1}, template::AbstractVector{T2}, element_type::Type{T3}=Float64; 
    normalize_template=true, direct=true) where {T1 <: Number, T2 <: Number, T3 <: AbstractFloat}
    
    if isempty(series)
        throw(ArgumentError("Series must be a non-empty vector."))
    elseif isempty(template)
        throw(ArgumentError("Template must be a non-empty vector."))
    elseif size(series) < size(template)
        throw(DimensionMismatch("Template is longer than series."))
    end

    if normalize_template
        template_mean, template_std = mean_and_std(template, corrected=false)
        template = @. (template - template_mean) / template_std
    end

    cc = similar(series, element_type, length(series) - length(template) + 1)
    if direct
        crosscorrelatedirect!(cc, series, template)
    else
        crosscorrelatefft!(cc, series, template)
    end
    cc
end

function crosscorrelatedirect!(cc, series, template)
    N = length(template)
    series_segment = similar(template) 
    for n in eachindex(cc)
        segment_view = view(series, n:n + N - 1)
        segment_std = std(segment_view, corrected=false)
        @. series_segment = segment_view / segment_std
        @inbounds cc[n] = dot(series_segment, template) / N
    end
end

function crosscorrelatefft!(cc, series, template)
    throw(ErrorException("Not yet implemented."))
end

"""
    maxfilter(x, tolerance)

Return a vector of the same size of `x` whose element at index `n` is the maximum 
of the elements of `x` which are at most `tolerance` apart from `x[n]` 
(i.e. in a window of at most `2 * tolerance + 1` elements).

# Examples

```jldoctest
julia> maxfilter(sin.(0:0.25pi:2pi), 1)
9-element Vector{Float64}:
  0.7071067811865475
  1.0
  1.0
  1.0
  0.7071067811865476
  1.2246467991473532e-16
 -0.7071067811865475
 -2.4492935982947064e-16
 -2.4492935982947064e-16
```
"""
function maxfilter(x::AbstractVector, l)
    y = similar(x)
    for n in eachindex(x)
        lower = max(n - l, firstindex(x))
        upper = min(n + l, lastindex(x))
        y[n] = maximum(view(x, lower:upper))
    end
    y
end

"""
    stack(correlations, offsets)

Return the average cross-correlation after aligning `correlations`. 
Each cross-correlation `correlations[n]` is shifted to the right by `offsets[n]`.

# Examples

```jldoctest
julia> stack([[0, 1.0, 0, 0], [0, 0, 1.0, 0]], [1, 2])
3-element OffsetArray(::Vector{Float64}, 0:2) with eltype Float64 with indices 0:2:
 0.0
 1.0
 0.0
```
"""
function stack(correlations::AbstractVector{T1}, offsets::AbstractVector{T2}) where {T1 <: AbstractVector, T2 <: Integer}
    if isempty(correlations)
        throw(ArgumentError("Cross-correlations vector must be non-empty."))
    elseif length(correlations) != length(offsets)
        throw(DimensionMismatch("Cross-correlations and offsets vectors must have the same length."))
    end
    stackedcorrelations = OffsetVector.(correlations, -offsets)
    start = maximum(series -> firstindex(series), stackedcorrelations)
    stop = minimum(series -> lastindex(series), stackedcorrelations)
    OffsetVector(mean(series -> view(series, start:stop), stackedcorrelations), start:stop)
end

"""
    correlatetemplate(data, template, offsets, tolerance, [element_type=Float64])

Return the cross-correlation between `data` and `template`. The cross-correlations for each series 
are aligned using `offsets` and averaged. If `tolerance`` is not zero, then the average accounts for possible
misplacement of each series by `tolerance` sample, and return the average of the maximum cross-correlation 
compatible with that misplacement.
"""
function correlatetemplate(data, template, offsets, tolerance, element_type=Float64)
    if isempty(data)
        throw(ArgumentError("Data must be non-empty."))
    elseif !(length(data) == length(template) == length(offsets))
        throw(DimensionMismatch("Data, template and shifts must have the same length."))
    end
    correlations = similar(data, Vector{element_type})
    Threads.@threads for n = eachindex(data)
        correlations[n] = maxfilter(crosscorrelate(data[n], template[n], element_type), tolerance)
    end
    stack(correlations, offsets)
end