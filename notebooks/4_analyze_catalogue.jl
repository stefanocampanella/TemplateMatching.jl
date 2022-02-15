### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ 95347531-11f3-4a66-ba6c-7ffdc70b7806
using JLD2

# ╔═╡ f07ab63f-dfb4-493c-a40c-6b3363a858c3
using Plots

# ╔═╡ 356cc4c2-e565-4847-9286-cfc1f835654d
using StatsBase

# ╔═╡ a7fcfdec-db16-4a40-b2f6-5c9c6f856ffc
using Printf

# ╔═╡ 9406c22e-fb52-461c-9b48-22bceb1b056c
md"## Reading templates catalogue"

# ╔═╡ 2fc4fdfe-81b2-45f1-9ed4-b8ad5bb86ee3
md"""Path of the CSV templates catalogue $(@bind templatespath TextField(default="../data/2021-01-12_20-25-30/cat2021-01-12_20-25-30_gab6_.2021.5.10.23.37.11.3267.csv"))"""

# ╔═╡ 6ac6cce3-09ae-421a-b502-b5c8807bf555
md"## Reading Sensors Positions"

# ╔═╡ 7d1f6c79-afbc-4300-bb6b-3aeed6844582
md"""Path of the JLD2 file containing continuous data: $(@bind datapath TextField(default="../data/2021-01-12_20-25-30/data.jld2"))"""

# ╔═╡ 05e33788-5c54-4425-b1ad-dba175ad5cbe
data = load(datapath, "data")

# ╔═╡ be090485-7917-464a-8e8c-aec4a2c009d3
begin
	sensors_positions = CSV.read("../data/2021-01-12_20-25-30/passive-xyz.csv", DataFrame, header=[:north, :east, :up])
	for s = [:north, :east, :up]
		sensors_positions[!, s] .*= 100
	end
	sensors_positions.sensor = axes(sensors_positions, 1) .- 1
	filter!(r -> r.sensor in keys(data), sensors_positions)
	sensors_positions
end

# ╔═╡ 57b35604-727b-4eac-bcbd-ea302e4c79a2
md"## Reading template matched catalogue"

# ╔═╡ b62d7318-3a6d-4442-b6cc-3166b87dff8a
md"""Path of the CSV template-match catalogue $(@bind cataloguepath TextField(default="../data/2021-01-12_20-25-30/templatematch.csv"))"""

# ╔═╡ c46997fc-9ab5-4678-bdd1-6874f77e7bca
catalogue = CSV.read(cataloguepath, DataFrame)

# ╔═╡ 102e4534-6be4-4436-be36-69726a393253
md" ## Statistics of the detections "

# ╔═╡ 45ef95ec-c638-4a94-b631-213cf3170fab
nontemplates = @. !(catalogue.correlation > 0.999 && catalogue.relative_magnitude ≈ 0.0)

# ╔═╡ feebcd5f-91c4-4a70-a064-e642e3172614
histogram(catalogue[nontemplates, :correlation], label=nothing)

# ╔═╡ 476101af-0108-4628-9ac1-f233a456a869
summarystats(catalogue[nontemplates, :correlation])

# ╔═╡ 3be237b0-9f19-4e9a-aa13-c6355771b8e8
histogram(catalogue[nontemplates, :relative_magnitude], label=nothing)

# ╔═╡ 7453767a-e8f5-4446-b198-fbe48c992aa5
summarystats(catalogue[nontemplates, :relative_magnitude])

# ╔═╡ 96cf6e3b-1c2b-49f6-979d-0fce28a55a06
starttime = DateTime(2021, 01, 12, 20, 25, 30)

# ╔═╡ 56c57ae9-5254-4dca-a793-bfdeee9a0ba5
begin
	columns = [:Year, :Month, :Day, :Hour, :Minute, :Second, :North, :East, :Up]
	df = CSV.read(templatespath, DataFrame, select=columns)
	sec = floor.(Int, df.Second)
	usec = round.(Int, 1e3 * (df.Second - sec))
	templates = DataFrame()
	templates.datetime = DateTime.(df.Year, df.Month, df.Day, df.Hour, df.Minute, sec, usec)
	templates.sample = map(x -> 1000 * x.value, @. (templates.datetime - starttime))
	templates.north = df.North
	templates.east = df.East
	templates.up = df.Up
	templates
end

# ╔═╡ 08f3bb2c-c80e-43ff-814b-e17ae5a912dc
scatter(catalogue[:, :sample],
		templates[catalogue.template, :sample],
		zcolor=catalogue.correlation,
		c=:heat,
		markersize=4exp.(catalogue.relative_magnitude),
		ylabel="Template origin [μs]",
		xlabel="Detection origin [μs]",
		colorbar=:right,
		colorbar_title="\nCross-correlation",
   		right_margin = 3Plots.mm,
		legend=nothing,
		dpi=500)

# ╔═╡ 1d576253-f719-4f21-8cc2-e794399a6ad5
let nframes = 500, size_r = 0.5, alpha_r = 1e-6, data_len = maximum(length, values(data))
	@gif for n = range(1, data_len, nframes)
		selection = catalogue[catalogue.sample .<= n, :]
		if !isempty(selection)
			count_upto = combine(groupby(selection, :template), nrow => :count)
			lastsample_upto = combine(groupby(selection, :template), :sample => maximum)
			alphas = @. exp10(alpha_r * (lastsample_upto[:, :sample_maximum] - n))
	    	scatter(templates[count_upto[!, :template], :east],
					templates[count_upto[!, :template], :north],
					templates[count_upto[!, :template], :up],
					markersize=size_r .* count_upto[:, :count],
					markeralpha=alphas,
					xlim=(0, 25),
					ylim=(0, 25),
					zlim=(0, 25),
					title=@sprintf("%.2f ms", 1e-3n),
					legend=nothing)
		end
	end
end

# ╔═╡ Cell order:
# ╠═aa854f02-46bf-4a2b-b1e9-5fb52b99925c
# ╠═8b4426b8-588b-49ad-83c5-a79ed698d704
# ╟─9406c22e-fb52-461c-9b48-22bceb1b056c
# ╠═3c07c33b-7f91-4e2d-a27a-09a8924f328c
# ╠═b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
# ╠═76e2e659-c708-4a08-8263-def6c3dfa930
# ╟─2fc4fdfe-81b2-45f1-9ed4-b8ad5bb86ee3
# ╠═56c57ae9-5254-4dca-a793-bfdeee9a0ba5
# ╟─6ac6cce3-09ae-421a-b502-b5c8807bf555
# ╠═95347531-11f3-4a66-ba6c-7ffdc70b7806
# ╟─7d1f6c79-afbc-4300-bb6b-3aeed6844582
# ╠═05e33788-5c54-4425-b1ad-dba175ad5cbe
# ╠═be090485-7917-464a-8e8c-aec4a2c009d3
# ╟─57b35604-727b-4eac-bcbd-ea302e4c79a2
# ╟─b62d7318-3a6d-4442-b6cc-3166b87dff8a
# ╠═c46997fc-9ab5-4678-bdd1-6874f77e7bca
# ╟─102e4534-6be4-4436-be36-69726a393253
# ╠═f07ab63f-dfb4-493c-a40c-6b3363a858c3
# ╠═356cc4c2-e565-4847-9286-cfc1f835654d
# ╠═45ef95ec-c638-4a94-b631-213cf3170fab
# ╟─feebcd5f-91c4-4a70-a064-e642e3172614
# ╟─476101af-0108-4628-9ac1-f233a456a869
# ╟─3be237b0-9f19-4e9a-aa13-c6355771b8e8
# ╟─7453767a-e8f5-4446-b198-fbe48c992aa5
# ╠═96cf6e3b-1c2b-49f6-979d-0fce28a55a06
# ╠═08f3bb2c-c80e-43ff-814b-e17ae5a912dc
# ╠═a7fcfdec-db16-4a40-b2f6-5c9c6f856ffc
# ╠═1d576253-f719-4f21-8cc2-e794399a6ad5
