### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ aa854f02-46bf-4a2b-b1e9-5fb52b99925c
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
end

# ╔═╡ 8b4426b8-588b-49ad-83c5-a79ed698d704
using PlutoUI

# ╔═╡ 3c07c33b-7f91-4e2d-a27a-09a8924f328c
using DataFrames

# ╔═╡ b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
using CSV

# ╔═╡ 76e2e659-c708-4a08-8263-def6c3dfa930
using Dates

# ╔═╡ 533d8e9e-40ff-4740-bbd9-8199cf0bcef6
using TemplateMatching

# ╔═╡ f07ab63f-dfb4-493c-a40c-6b3363a858c3
using Plots

# ╔═╡ 356cc4c2-e565-4847-9286-cfc1f835654d
using StatsBase

# ╔═╡ a7fcfdec-db16-4a40-b2f6-5c9c6f856ffc
# ╠═╡ disabled = true
#=╠═╡
using Printf
  ╠═╡ =#

# ╔═╡ 33ff9417-b0b4-4f60-8c81-22c1181b37aa
using MultivariateStats

# ╔═╡ 3ec5251a-2e3e-4344-bb1f-fcd79daf9714
using LinearAlgebra

# ╔═╡ 9406c22e-fb52-461c-9b48-22bceb1b056c
md"## Reading templates catalogue"

# ╔═╡ 2fc4fdfe-81b2-45f1-9ed4-b8ad5bb86ee3
md"""Path of the CSV templates catalogue $(@bind templatespath TextField(default="../data/2021-01-12_20-25-30/cat.csv"))"""

# ╔═╡ 60ad8eb3-a313-4dce-82bd-a74a4d6842fd
begin 
	import Base.+, Base.show
	
	+(x::Dates.UTInstant{Microsecond}, d::Microsecond) = Dates.UTInstant(x.periods + d)

	function show(io::IO, ::MIME"text/plain", x::Dates.UTInstant{Microsecond})
		date = DateTime(Dates.UTM(round(Int, 1e-3 * x.periods.value)))
		show(io, MIME"text/plain"(), date)
	end
end

# ╔═╡ 74063783-3454-4e4b-95d5-0427edbe9380
function DateTimeMicrosecond(y, m, d, h, mi, s, us)
	rata = us + 1_000_000 * (s + 60mi + 3600h + 86400 * Dates.totaldays(y, m, d))
	Dates.UTInstant{Microsecond}(Microsecond(rata))
end

# ╔═╡ 171de4ad-485a-4f92-818a-50cfbfa2330b
starttime = DateTimeMicrosecond(2021, 01, 12, 20, 25, 30, 0)

# ╔═╡ c76f6541-b0e2-4aba-9533-8ab6c23685bb
samplefreq = 2 # MHz

# ╔═╡ 56c57ae9-5254-4dca-a793-bfdeee9a0ba5
templates = let
	columns = [:Year, :Month, :Day, :Hour, :Minute, :Second, :North, :East, :Up]
	df = CSV.read(templatespath, DataFrame, select=columns)
	sec = floor.(Int, df.Second)
	usec = round.(Int, 1e6 * (df.Second - sec))
	templates = DataFrame()
	templates.datetime = DateTimeMicrosecond.(df.Year, df.Month, df.Day, df.Hour, df.Minute, sec, usec)
	templates.origin_sample = map(x -> round(Int, x.value * samplefreq), templates.datetime .- starttime)
	templates.origin_time = templates.origin_sample ./ samplefreq
	templates.north = df.North
	templates.east = df.East
	templates.up = df.Up
	templates
end

# ╔═╡ 57b35604-727b-4eac-bcbd-ea302e4c79a2
md"## Reading template matched catalogue"

# ╔═╡ b62d7318-3a6d-4442-b6cc-3166b87dff8a
md"""Path of the CSV template-match catalogue $(@bind cataloguepath TextField(default="../data/2021-01-12_20-25-30/augmented_catalogue.csv"))"""

# ╔═╡ c46997fc-9ab5-4678-bdd1-6874f77e7bca
catalogue = CSV.read(cataloguepath, DataFrame)

# ╔═╡ 08f3bb2c-c80e-43ff-814b-e17ae5a912dc
scatter(catalogue[:, :origin_time],
		templates[catalogue.template, :origin_time],
		zcolor=catalogue.crosscorrelation,
		c=:heat,
		markersize=4exp.(catalogue.relative_magnitude),
		ylabel="Template origin [μs]",
		xlabel="Detection origin [μs]",
		colorbar=:right,
		colorbar_title="\nCross-correlation",
   		right_margin = 3Plots.mm,
		legend=nothing,
		dpi=500)

# ╔═╡ 102e4534-6be4-4436-be36-69726a393253
md" ## Statistics of the detections "

# ╔═╡ 2ce75a18-b5d6-4e6f-81ce-875527b14092
selection = let
	# Skip matches without enough channels
	relocatable_matches = dropmissing(catalogue)
	# Skip matches with residuals too large or small
	relocatable_matches[.!TemplateMatching.mad_test(relocatable_matches.multilateration_residual), :]
end

# ╔═╡ 567cbd00-b375-417f-a7f5-8ac9a4bc4c7f
nontemplates = @. !(selection.crosscorrelation > 0.999 && selection.relative_magnitude ≈ 0.0)

# ╔═╡ feebcd5f-91c4-4a70-a064-e642e3172614
histogram(selection[nontemplates, :crosscorrelation], label=nothing)

# ╔═╡ 476101af-0108-4628-9ac1-f233a456a869
summarystats(selection[nontemplates, :crosscorrelation])

# ╔═╡ 3be237b0-9f19-4e9a-aa13-c6355771b8e8
histogram(selection[nontemplates, :relative_magnitude], label=nothing)

# ╔═╡ 7453767a-e8f5-4446-b198-fbe48c992aa5
summarystats(selection.relative_magnitude)

# ╔═╡ 2f98cf67-87ee-4f12-a507-b2259d780ffd
md"## Animations"

# ╔═╡ 6dddcebd-d730-470e-a831-154ee3102e0d
#=╠═╡
let nframes = 500, size_r = 0.8, alpha_r = 1e-6, data_len = Int(4e6)
	@gif for n = range(1, data_len, nframes)
		selection_upto = selection[selection.origin_time .<= n / samplefreq, :]
		if !isempty(selection_upto)
	    	scatter(selection_upto[!, :east],
					selection_upto[!, :north],
					selection_upto[!, :up],
					markersize=2.5,
					xlim=(0, 25),
					ylim=(0, 25),
					zlim=(0, 25),
					title=@sprintf("%.2f ms", 1e-3n / samplefreq),
					legend=nothing)
		end
	end
end
  ╠═╡ =#

# ╔═╡ 1d576253-f719-4f21-8cc2-e794399a6ad5
# ╠═╡ disabled = true
#=╠═╡
let nframes = 500, size_r = 0.8, alpha_r = 1e-6, data_len = Int(4e6)
	@gif for n = range(1, data_len, nframes)
		selection_upto = selection[selection.sample .<= n, :]
		if !isempty(selection_upto)
			count_upto = combine(groupby(selection_upto, :template), nrow => :count)
			lastsample_upto = combine(groupby(selection_upto, :template), :sample => maximum)
			alphas = @. exp10(alpha_r * (lastsample_upto[:, :sample_maximum] - n))
	    	scatter(templates[count_upto[!, :template], :east],
					templates[count_upto[!, :template], :north],
					templates[count_upto[!, :template], :up],
					markersize=size_r .* count_upto[:, :count],
					markeralpha=alphas,
					xlim=(0, 25),
					ylim=(0, 25),
					zlim=(0, 25),
					title=@sprintf("%.2f ms", 1e-3n / samplefreq),
					legend=nothing)
		end
	end
end
  ╠═╡ =#

# ╔═╡ 0d237fe0-35bb-4624-a49a-6086c5d9ed09
md"## Principal Component Analysis"

# ╔═╡ 0131b25b-8b5c-4d1e-a520-830f5f20f345
X = transpose(Matrix(selection[!, [:north, :east, :up]]))

# ╔═╡ 09731759-c6c6-43ba-bc99-493327efa342
M = fit(PCA, X, maxoutdim=2)

# ╔═╡ 8c740959-ee9b-4139-a31b-fdd7ccc7f563
xlim, ylim = extrema(predict(M, X), dims=2)

# ╔═╡ 2837928e-692a-42ee-a064-07a3042d3d7b
# ╠═╡ disabled = true
#=╠═╡
let nframes = 500, data_len = Int(4e6)
	@gif for n = range(1, data_len, nframes)
		selection_upto = selection[selection.origin_time .<= n / samplefreq, :]
		if !isempty(selection)
			X = transpose(Matrix(selection_upto[!, [:north, :east, :up]]))
			Y = predict(M, X)
	    	scatter(Y[1, :],
					Y[2, :],
					xlim=xlim,
					ylim=ylim,
					title=@sprintf("%.2f ms", 1e-3n / samplefreq),
					legend=nothing)			
		end
	end
end
  ╠═╡ =#

# ╔═╡ 3212d37b-e498-4711-a067-14dca93a91a8
let
	Y = predict(M, X)
	plot(fit(Histogram, (Y[1, :], Y[2, :]), nbins=64))
end

# ╔═╡ 48bf7fda-607b-44fe-85aa-4e57ad9c565a
projection(M)

# ╔═╡ 04d6e989-3fd0-4efd-a570-0ae34164cf4f
let
	Y = predict(M, X)
	Z = norm.(eachcol(X .- reconstruct(M, Y)))
	selection = Z .< percentile(Z, 90)
	scatter(Y[1, selection], Y[2, selection], zcolor=Z[selection], c=:heat, markersize=2.5, label="")
end

# ╔═╡ Cell order:
# ╠═aa854f02-46bf-4a2b-b1e9-5fb52b99925c
# ╠═8b4426b8-588b-49ad-83c5-a79ed698d704
# ╟─9406c22e-fb52-461c-9b48-22bceb1b056c
# ╠═3c07c33b-7f91-4e2d-a27a-09a8924f328c
# ╠═b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
# ╠═76e2e659-c708-4a08-8263-def6c3dfa930
# ╠═2fc4fdfe-81b2-45f1-9ed4-b8ad5bb86ee3
# ╠═74063783-3454-4e4b-95d5-0427edbe9380
# ╠═171de4ad-485a-4f92-818a-50cfbfa2330b
# ╠═60ad8eb3-a313-4dce-82bd-a74a4d6842fd
# ╠═c76f6541-b0e2-4aba-9533-8ab6c23685bb
# ╠═56c57ae9-5254-4dca-a793-bfdeee9a0ba5
# ╟─57b35604-727b-4eac-bcbd-ea302e4c79a2
# ╠═b62d7318-3a6d-4442-b6cc-3166b87dff8a
# ╠═533d8e9e-40ff-4740-bbd9-8199cf0bcef6
# ╠═f07ab63f-dfb4-493c-a40c-6b3363a858c3
# ╠═c46997fc-9ab5-4678-bdd1-6874f77e7bca
# ╠═08f3bb2c-c80e-43ff-814b-e17ae5a912dc
# ╟─102e4534-6be4-4436-be36-69726a393253
# ╠═356cc4c2-e565-4847-9286-cfc1f835654d
# ╠═2ce75a18-b5d6-4e6f-81ce-875527b14092
# ╠═567cbd00-b375-417f-a7f5-8ac9a4bc4c7f
# ╠═feebcd5f-91c4-4a70-a064-e642e3172614
# ╠═476101af-0108-4628-9ac1-f233a456a869
# ╠═3be237b0-9f19-4e9a-aa13-c6355771b8e8
# ╟─7453767a-e8f5-4446-b198-fbe48c992aa5
# ╟─2f98cf67-87ee-4f12-a507-b2259d780ffd
# ╠═a7fcfdec-db16-4a40-b2f6-5c9c6f856ffc
# ╠═6dddcebd-d730-470e-a831-154ee3102e0d
# ╠═1d576253-f719-4f21-8cc2-e794399a6ad5
# ╟─0d237fe0-35bb-4624-a49a-6086c5d9ed09
# ╠═33ff9417-b0b4-4f60-8c81-22c1181b37aa
# ╠═0131b25b-8b5c-4d1e-a520-830f5f20f345
# ╠═09731759-c6c6-43ba-bc99-493327efa342
# ╠═8c740959-ee9b-4139-a31b-fdd7ccc7f563
# ╠═2837928e-692a-42ee-a064-07a3042d3d7b
# ╠═3ec5251a-2e3e-4344-bb1f-fcd79daf9714
# ╠═3212d37b-e498-4711-a067-14dca93a91a8
# ╠═48bf7fda-607b-44fe-85aa-4e57ad9c565a
# ╠═04d6e989-3fd0-4efd-a570-0ae34164cf4f
