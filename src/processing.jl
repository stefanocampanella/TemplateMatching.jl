
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
                push!(peaks, midpoint)
                push!(heights, height)
            end
        end
    end

    if threshold > 0 || distance > 0
        tokeep = trues(length(peaks))
        if threshold > 0
            tokeep .&= selectbypeaksheight(heights, threshold)
        end
        if distance > 0
            tokeep .&= selectbypeaksdistance(peaks, heights, distance)
        end
        peaks[tokeep], heights[tokeep]
    else
        peaks, heights
    end
end

function selectbypeaksheight(heights, threshold)
    heights .> threshold
end

function selectbypeaksdistance(peaks, heights, distance)
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

    tokeep
end

"""
    mad_test(xs, [r=3])

Return a vector of booleans whose elements are `true` if the absolute deviation 
from the median of `xs` of the corresponding elements (of `xs`) is at least 
`r` times the median absolute deviation, and `false` otherwise.

# Examples
```jldoctest
julia> TemplateMatching.mad_test([1, 2, 3, 100])
4-element BitVector:
 0
 0
 0
 1
```
"""
function mad_test(xs, r=3.0)
    deviations = abs.(xs .- median(xs))
    deviations .> r * median(deviations)
end

function series_magnitude(series, template, index)
    series_amp = maximum(abs.(view(series, index:index + size(template, 1) - 1)))
    template_amp = maximum(abs.(template))
    if iszero(template_amp)
        NaN
    elseif iszero(series_amp)
        NaN
    else
        log10(series_amp / template_amp)
    end
end

"""
    magnitude(data, template, indices, [outlier=TemplateMatching.mad], [avg=StatsBase.mean])

Return the average relative magnitude of the view of series in `data` starting at `indices` and 
having the same length of the corresponding `template` series. The outliers, i.e. the elements 
corresponding to those of `outlier(magnitudes)` evaluates that are true, are excluded. The average
is calculated using `avg(magnitudes)`.
"""
function magnitude(data, template, indices, outlier=mad_test, avg=mean)
    magnitudes = series_magnitude.(data, template, indices)
    avg(magnitudes[.!outlier(magnitudes)])
end

"""
    findmax_window(crosscorrelation, index, tolerance)

Return a tuples containing the maximum of `crosscorrelation`, in a window 
centered around `index` and radius equal to `tolerace`, and its index. 

# Examples

```jldoctest
julia> findmax_window.([abs.(sin.(0:0.01pi:4pi)) for _ = 1:4], [49, 153, 249, 353], 5)
4-element Vector{Tuple{Float64, Int64}}:
 (1.0, 51)
 (1.0, 151)
 (1.0, 251)
 (1.0, 351)
```
"""
function findmax_window(crosscorrelation, index, tolerance)
    if isempty(crosscorrelation)
        throw(ArgumentError("Cross-correlation series must be non-empty."))
    elseif !(index in axes(crosscorrelation, 1))
        throw(ArgumentError("Index must be a valid index of cross-correlation series."))
    end
    start = max(index  - tolerance, first(axes(crosscorrelation, 1)))
    stop = min(index  + tolerance, last(axes(crosscorrelation, 1)))
    cc_max, indx_max = findmax(view(crosscorrelation, start:stop))
    cc_max, start + indx_max - 1
end
