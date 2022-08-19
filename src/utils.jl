"""
    subsampleshift(waveform, quantity, [algorithm=BSpline(Cubic(Interpolations.Flat(OnGrid())))])

Return `waveform` shifted by `quantity` using `algorithm` for interpolation.

# Examples

```jldoctest
julia> let f(x) = x^3 + 2x^2 + 3x + 4, x = range(0., 1., 256)
       isapprox(f.(x), subsampleshift(subsampleshift(f.(x), 0.5), -0.5), atol=1e-2)
       end
true
```
"""
function subsampleshift(waveform::AbstractVector{<:AbstractFloat}, quantity; algorithm=BSpline(Cubic(Interpolations.Flat(OnGrid()))))
    if abs(quantity) > 1.0
        throw(ArgumentError("Shift should be less than 1.0, got $quantity"))
    else
        xs = axes(waveform, 1)
        itp = interpolate(waveform, algorithm)
        ext = extrapolate(itp, Interpolations.Flat())
        @. ext(xs + quantity)
    end
end


"""
    findmax_window(crosscorrelation, index, tolerance)

Return a tuple containing the maximum of `crosscorrelation`, in a window 
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
    start = max(index  - tolerance, firstindex(crosscorrelation))
    stop = min(index  + tolerance, lastindex(crosscorrelation))
    if start == firstindex(crosscorrelation) && stop == lastindex(crosscorrelation)
        findmax(crosscorrelation)
    else
        cc_max, indx_max = findmax(view(crosscorrelation, start:stop))
        cc_max, start + indx_max - 1
    end
end


"""
    estimatetoa(trace, waveform, center, left, tolerance)

Return a tuple containing the estimated position, with subsample precision, 
and value of maximum cross-correlation between `trace` and `waveform` 
in a window around `center` and radius equal to `tolerace`. 
The estimation is performed by approximating the cross-correlation with
a parabola near the maximum. This estimate is biased and its accuracy 
decreases with decreasing sampling frequencies.
The position of the maximum cross-correlation is the Time Of Arrival (TOA) 
in sample units, i.e. the true TOA divided by sampling frequency.

# Examples

```jldoctest
julia> let v = sin.(range(0, pi, 256)), x = [zeros(128); v; zeros(128)], y = TemplateMatching.subsampleshift(v, 0.5)
       toa, cc = estimatetoa(x, y, 128, 5)
       round(toa, digits=1), round(cc, digits=1)
       end
(129.5, 1.0)
```
"""
function estimatetoa(trace, waveform, center, tol)
    a, b = max(firstindex(trace), center - tol), min(lastindex(trace), center + tol + length(waveform) - 1)
    cc = OffsetVector(crosscorrelate(view(trace, a:b), waveform, usefft=false), center - tol - 1)
    cc_max, n_max = findmax_window(cc, center, tol)
    if firstindex(cc) < n_max < lastindex(cc) && cc_max < 1.0 && cc_max ≉ 1.0
        y1, y2, y3 = cc[n_max - 1:n_max + 1]
        delta = (y1 - y3) / 2(y1 - 2y2 + y3)
        n_max + delta, min(1.0, y2 + (y3 - y1) / 4delta)
    else
        Float64(n_max), cc_max
    end
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
function mad_test(xs::AbstractVector{<:Number}, r::Number)
    if isempty(xs)
        weights(BitVector())
    else
        deviations = abs.(xs .- median(xs))
        weights(deviations .> r * median(deviations))
    end
end


"""
    relative_magnitude(trace_a, trace_b)

Return the log10 of the ratio of peak amplitudes of `trace_a` and `trace_b`.
"""
function relative_magnitude(trace_a::AbstractVector{T}, trace_b::AbstractVector{T}) where T <: Number
    if isempty(trace_a) || isempty(trace_b)
        throw(ArgumentError("Arguments must be non empty vectors."))
    elseif length(trace_a) != length(trace_b)
        throw(ArgumentError("Arguments must be vectors of the same length."))

    else
        trace_a_peak_amplitude = maximum(abs.(trace_a))
        trace_b_peak_amplitude = maximum(abs.(trace_b))
        if iszero(trace_b_peak_amplitude)
            NaN
        elseif iszero(trace_a_peak_amplitude)
            NaN
        else
            log10(trace_a_peak_amplitude / trace_b_peak_amplitude)
        end
    end
end