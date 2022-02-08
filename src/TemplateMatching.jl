module TemplateMatching

using StatsBase
using LinearAlgebra
using OffsetArrays

export crosscorrelate, maxfilter, stack, correlatetemplate, findpeaks

"""
    crosscorrelate(series, template, [cc_eltype=Float64], normalize_template=true)

Return the normalized cross-correlation between `series` and `template`. 

If `series` is a vector of ``n`` elements and template a vector of ``N`` elements, 
cross-correlation will be a vector of ``n - N + 1`` elements of type `cc_eltype`. 
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
function crosscorrelate(series::AbstractVector{T1}, template::AbstractVector{T2}, cc_eltype::Type{T3}=Float64; 
    normalize_template=true) where {T1 <: Number, T2 <: Number, T3 <: AbstractFloat}
    
    if isempty(series)
        throw(ArgumentError("Series must be a non-empty vector."))
    elseif isempty(template)
        throw(ArgumentError("Template must be a non-empty vector."))
    elseif size(series, 1) < size(template, 1)
        throw(DimensionMismatch("Template is longer than series."))
    end

    if normalize_template
        template_mean, template_std = mean_and_std(template, corrected=false)
        y = @. ($collect(template) - template_mean) / template_std
    else
        y = collect(template)
    end

    N = size(template, 1)
    cc = Vector{cc_eltype}(undef, size(series, 1) - N + 1)
    x = similar(y) 
    for n in eachindex(cc)
        series_view = view(collect(series), n:n + N - 1)
        series_view_std = std(series_view, corrected=false)
        @. x = series_view / series_view_std
        cc[n] = dot(x, y) / N
    end
    cc
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
function maxfilter(x, l)
    y = similar(x)
    for n in eachindex(x)
        lower = max(n - l, first(axes(x, 1)))
        upper = min(n + l, last(axes(x, 1)))
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
    elseif size(correlations, 1) != size(offsets, 1)
        throw(DimensionMismatch("Cross-correlations and offsets vectors must have the same length."))
    end
    stackedcorrelations = OffsetVector.(correlations, -offsets)
    start = maximum(series -> first(axes(series, 1)), stackedcorrelations)
    stop = minimum(series -> last(axes(series, 1)), stackedcorrelations)
    OffsetVector(mean(series -> view(series, start:stop), stackedcorrelations), start:stop)
end

"""
    correlatetemplate(data, template, offsets, tolerance, [cc_eltype=Float64])

Return the cross-correlation between `data` and `template`. The cross-correlations for each series 
are aligned using `offsets` and averaged. If `tolerance`` is not zero, then the average accounts for possible
misplacement of each series by `tolerance` sample, and return the average of the maximum cross-correlation 
compatible with that misplacement.
"""
function correlatetemplate(data, template, offsets, tolerance, cc_eltype=Float64)
    if isempty(data)
        throw(ArgumentError("Data must be non-empty."))
    elseif !(size(data, 1) == size(template, 1) == size(offsets, 1))
        throw(DimensionMismatch("Data, template and shifts must have the same length."))
    end
    correlations = similar(data, Vector{cc_eltype})
    Threads.@threads for n = eachindex(data)
        correlations[n] = maxfilter(crosscorrelate(data[n], template[n], cc_eltype), tolerance)
    end
    stack(correlations, offsets)
end

"""
    findpeaks(signal, threshold, distance)
Find all the local maxima of `signal` above `threshold` and at least `distance` apart. 
The maxima closer than `distance` are removed starting from the ones with 
higher neighbour maxima up until the condition is met. 

Return two vectors one with indices and one with heights of the peaks.

The implementation is taken from [`scipy.signal`](https://github.com/scipy/scipy/blob/\
8a64c938ddf1ae4c02a08d2c5e38daeb8d061d38/scipy/signal/_peak_finding.py#L723-L1003).

# Examples

```jldoctest
julia> findpeaks(abs.(sin.(0:0.01pi:2pi)), 0.5, 0)
([51, 151], [1.0, 1.0])
```
"""
function findpeaks(signal::AbstractVector, threshold::Number, distance::Number)
	start = first(axes(signal, 1)) + 1
	stop = last(axes(signal, 1)) - 1
	peaks = Vector{eltype(axes(signal, 1))}(undef, 0)
	heights = Vector{eltype(signal)}(undef, 0)
	for n in start:stop
		if signal[n-1] < signal[n]
			height = signal[n]
			n_ahead = n + 1
			while height == signal[n_ahead] && n_ahead < stop
				n_ahead += 1
			end
			if signal[n_ahead] < height
				midpoint = div(n + n_ahead, 2)
				if height > threshold
					push!(peaks, midpoint)
					push!(heights, height)
				end
			end
		end
	end
	
	if distance > 0
		npeaks = length(peaks)
		tokeep = fill(true, npeaks)
		for j in sortperm(heights, rev=true)

			if !tokeep[j]
				continue
			end

			k = j - 1
			while k >= 1 && (peaks[j] - peaks[k]) <= distance
				tokeep[k] = false
				k -= 1
			end

			k = j + 1
			while k <= npeaks && (peaks[k] - peaks[j]) <= distance
				tokeep[k] = false
				k += 1
			end
		end
	
		peaks[tokeep], heights[tokeep]
	else
		peaks, heights
	end
end

end
