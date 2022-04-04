var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TemplateMatching","category":"page"},{"location":"#TemplateMatching","page":"Home","title":"TemplateMatching","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TemplateMatching.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TemplateMatching]","category":"page"},{"location":"#TemplateMatching.correlatetemplate","page":"Home","title":"TemplateMatching.correlatetemplate","text":"correlatetemplate(data, template, offsets, tolerance, [element_type=Float64])\n\nReturn the cross-correlation between data and template. The cross-correlations for each series  are aligned using offsets and averaged. If toleranceis not zero, then the average accounts for possible misplacement of each series bytolerance` sample, and return the average of the maximum cross-correlation  compatible with that misplacement.\n\n\n\n\n\n","category":"function"},{"location":"#TemplateMatching.crosscorrelate-Union{Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{AbstractVector{T1}, AbstractVector{T2}}, Tuple{AbstractVector{T1}, AbstractVector{T2}, Type{T3}}} where {T1<:Number, T2<:Number, T3<:AbstractFloat}","page":"Home","title":"TemplateMatching.crosscorrelate","text":"crosscorrelate(series, template, [element_type=Float64], normalize_template=true)\n\nReturn the normalized cross-correlation between series and template. \n\nIf series is a vector of n elements and template a vector of N elements,  cross-correlation will be a vector of n - N + 1 elements of type element_type.  The element of cross-correlation χ_k is Pearson correlation coefficient between  x^(k)_i denoting the elements of series x_i + k and  the elements of template y_i, with i = 1 2 dots N\n\nχ_k = fracsumleft(x^(k)_i - μ_x^(k)right)  left(y_i - μ_yright)N σ_x^(k) σ_y  \n\nwhere μ and σ denotes the mean and variance respectively. \n\nIf normalize_template is set to false then template is assumed to have mean and std respectively equal to 0.0 and 1.0.\n\nSee the Wikipedia page  on cross-correlation for more details.\n\nExamples\n\njulia> crosscorrelate(sin.(0:0.25pi:2pi), [1, 1+√2, 1])\n7-element Vector{Float64}:\n  0.23258781949447394\n  1.000000000000001\n  0.23258781949447402\n  7.401486830834377e-17\n -0.23258781949447394\n -1.0000000000000007\n -0.23258781949447394\n\n\n\n\n\n","category":"method"},{"location":"#TemplateMatching.estimatetoa-NTuple{4, Any}","page":"Home","title":"TemplateMatching.estimatetoa","text":"estimatetoa(trace, waveform, center, left, tolerance)\n\nReturn a tuple containing the estimated position, with subsample precision,  and value of maximum cross-correlation between trace and waveform  in a window around center and radius equal to tolerace.  The estimation is performed by approximating the cross-correlation with a parabola near the maximum. This estimate is biased and its accuracy  decreases with decreasing sampling frequencies. The position of the maximum cross-correlation is the Time Of Arrival (TOA)  in sample units, i.e. the true TOA divided by sampling frequency.\n\n\n\n\n\n","category":"method"},{"location":"#TemplateMatching.findmax_window-Tuple{Any, Any, Any}","page":"Home","title":"TemplateMatching.findmax_window","text":"findmax_window(crosscorrelation, index, tolerance)\n\nReturn a tuple containing the maximum of crosscorrelation, in a window  centered around index and radius equal to tolerace, and its index. \n\nExamples\n\njulia> findmax_window.([abs.(sin.(0:0.01pi:4pi)) for _ = 1:4], [49, 153, 249, 353], 5)\n4-element Vector{Tuple{Float64, Int64}}:\n (1.0, 51)\n (1.0, 151)\n (1.0, 251)\n (1.0, 351)\n\n\n\n\n\n","category":"method"},{"location":"#TemplateMatching.findpeaks-Tuple{AbstractVector, Number, Number}","page":"Home","title":"TemplateMatching.findpeaks","text":"findpeaks(signal, threshold, distance)\n\nFind all the local maxima of signal above threshold and at least distance apart.  The maxima closer than distance are removed starting from the ones with  higher neighbour maxima up until the condition is met. \n\nReturn two vectors one with indices and one with heights of the peaks.\n\nThe implementation is taken from scipy.signal.\n\nExamples\n\njulia> findpeaks(abs.(sin.(0:0.01pi:2pi)), 0.5, 0)\n([51, 151], [1.0, 1.0])\n\n\n\n\n\n","category":"method"},{"location":"#TemplateMatching.mad_test","page":"Home","title":"TemplateMatching.mad_test","text":"mad_test(xs, [r=3])\n\nReturn a vector of booleans whose elements are true if the absolute deviation  from the median of xs of the corresponding elements (of xs) is at least  r times the median absolute deviation, and false otherwise.\n\nExamples\n\njulia> TemplateMatching.mad_test([1, 2, 3, 100])\n4-element BitVector:\n 0\n 0\n 0\n 1\n\n\n\n\n\n","category":"function"},{"location":"#TemplateMatching.magnitude","page":"Home","title":"TemplateMatching.magnitude","text":"magnitude(data, template, indices, [outlier=TemplateMatching.mad], [avg=StatsBase.mean])\n\nReturn the average relative magnitude of the view of series in data starting at indices and  having the same length of the corresponding template series. The outliers, i.e. the elements  corresponding to those of outlier(magnitudes) evaluates that are true, are excluded. The average is calculated using avg(magnitudes).\n\n\n\n\n\n","category":"function"},{"location":"#TemplateMatching.maxfilter-Tuple{AbstractVector, Any}","page":"Home","title":"TemplateMatching.maxfilter","text":"maxfilter(x, tolerance)\n\nReturn a vector of the same size of x whose element at index n is the maximum  of the elements of x which are at most tolerance apart from x[n]  (i.e. in a window of at most 2 * tolerance + 1 elements).\n\nExamples\n\njulia> maxfilter(sin.(0:0.25pi:2pi), 1)\n9-element Vector{Float64}:\n  0.7071067811865475\n  1.0\n  1.0\n  1.0\n  0.7071067811865476\n  1.2246467991473532e-16\n -0.7071067811865475\n -2.4492935982947064e-16\n -2.4492935982947064e-16\n\n\n\n\n\n","category":"method"},{"location":"#TemplateMatching.stack-Union{Tuple{T2}, Tuple{T1}, Tuple{AbstractVector{T1}, AbstractVector{T2}}} where {T1<:(AbstractVector), T2<:Integer}","page":"Home","title":"TemplateMatching.stack","text":"stack(correlations, offsets)\n\nReturn the average cross-correlation after aligning correlations.  Each cross-correlation correlations[n] is shifted to the right by offsets[n].\n\nExamples\n\njulia> stack([[0, 1.0, 0, 0], [0, 0, 1.0, 0]], [1, 2])\n3-element OffsetArray(::Vector{Float64}, 0:2) with eltype Float64 with indices 0:2:\n 0.0\n 1.0\n 0.0\n\n\n\n\n\n","category":"method"}]
}
