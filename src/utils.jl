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
    cc = OffsetVector(crosscorrelate(view(trace, a:b), waveform), center - tol - 1)
    cc_max, n_max = findmax_window(cc, center, tol)
    if firstindex(cc) < n_max < lastindex(cc) && cc_max < 1.0 && cc_max ≉ 1.0
        y1, y2, y3 = cc[n_max - 1:n_max + 1]
        delta = (y1 - y3) / 2(y1 - 2y2 + y3)
        n_max + delta, min(1.0, y2 + (y3 - y1) / 4delta)
    else
        Float64(n_max), cc_max
    end
end

function line_element(xt, v)
    x = view(xt, 1:3)
    t = xt[4]
    dot(x, x) - v^2 * t^2
end

residue_rms(xt, sensors_readings_itr, v) = sqrt(mean(ys -> line_element(ys - xt, v)^2, sensors_readings_itr))

locate(sensors_readings_itr, v, guess) = optimize(xt -> residue_rms(xt, sensors_readings_itr, v), guess)
