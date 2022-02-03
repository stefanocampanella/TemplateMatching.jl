var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TemplateMatching","category":"page"},{"location":"#TemplateMatching","page":"Home","title":"TemplateMatching","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TemplateMatching.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TemplateMatching]","category":"page"},{"location":"#TemplateMatching.crosscorrelate-Union{Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{AbstractVector{T1}, AbstractVector{T2}}, Tuple{AbstractVector{T1}, AbstractVector{T2}, Type{T3}}} where {T1<:Number, T2<:Number, T3<:AbstractFloat}","page":"Home","title":"TemplateMatching.crosscorrelate","text":"crosscorrelate(series, template, cc_eltype=Float64; normalize_template=false) -> Vector{cc_eltype}\n\nReturn the normalized cross-correlation `cc` between `series` and `template`. If `series` is a vector of ``n`` elements \nand template a vector of ``k`` elements, then `cc` is a vector of ``n - k + 1`` elements of type `cc_eltype`. \nThe element of `cc` at position ``l`` is Pearson correlation coefficient between the segment of `series` \nstarting at ``l`` and ``k`` elements long. See the [Wikipedia page](https://en.wikipedia.org/wiki/Cross-correlation#Normalization)\non cross-correlation for more details.\n\n\n\n\n\n","category":"method"}]
}
