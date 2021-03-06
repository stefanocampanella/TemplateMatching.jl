
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
    start = firstindex(signal) + 1
    stop = lastindex(signal) - 1
    peaks = Vector{eltype(axes(signal, 1))}(undef, 0)
    heights = Vector{eltype(signal)}(undef, 0)
    for n in start:stop
        if signal[n-1] < signal[n]
            height = signal[n]
            n_ahead = n + 1
            while height == signal[n_ahead] && n_ahead < stop
                n_ahead += 1
            end
            if signal[n_ahead] < height && height > threshold
                midpoint = div(n + n_ahead, 2)
                push!(peaks, midpoint)
                push!(heights, height)
            end
        end
    end

    if distance > 0
        tokeep = selectbypeaksdistance(peaks, heights, distance)
        peaks[tokeep], heights[tokeep]
    else
        peaks, heights
    end
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
function mad_test(xs::AbstractVector{<:Number}, r=3.0)
    if isempty(xs)
        BitVector()
    else
        deviations = abs.(xs .- median(xs))
        deviations .> r * median(deviations)
    end
end

function series_magnitude(series, waveform, index)
    series_amp = maximum(abs.(view(series, index:index + length(waveform) - 1)))
    template_amp = maximum(abs.(waveform))
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
