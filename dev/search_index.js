var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TemplateMatching","category":"page"},{"location":"#TemplateMatching","page":"Home","title":"TemplateMatching","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TemplateMatching.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TemplateMatching]","category":"page"},{"location":"#TemplateMatching.crosscorrelate-Union{Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{AbstractVector{T1}, AbstractVector{T2}}, Tuple{AbstractVector{T1}, AbstractVector{T2}, Type{T3}}} where {T1<:Number, T2<:Number, T3<:AbstractFloat}","page":"Home","title":"TemplateMatching.crosscorrelate","text":"crosscorrelate(series, template, cc_eltype=Float64; normalize_template=false)\n\nReturn the normalized cross-correlation between series and template. \n\nIf series is a vector of n elements and template a vector of N elements,  cross-correlation will be a vector of n - N + 1 elements of type cc_eltype.  The element of cross-correlation χ_k is Pearson correlation coefficient between  x^k_i denoting the elements of series x_i + k and  the elements of template y_i, with i = 1 2 dots N\n\nχ_k = fracsumleft(x^k_i - μ_x^kright) left(y_i - μ_yright)N σ_x^k σ_y  \n\nwhere μ and σ denotes the mean and variance respectively.\n\nSee the Wikipedia page  on cross-correlation for more details.\n\nExamples\n\njulia> crosscorrelate(sin.(0:0.25pi:2pi), [1, 1+√2, 1])\n7-element Vector{Float64}:\n  0.23258781949447394\n  1.000000000000001\n  0.23258781949447402\n  7.401486830834377e-17\n -0.23258781949447394\n -1.0000000000000007\n -0.23258781949447394\n\n\n\n\n\n","category":"method"}]
}
