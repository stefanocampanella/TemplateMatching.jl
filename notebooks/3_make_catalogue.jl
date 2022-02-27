### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ d47c1cf4-4844-4c66-8a21-b62ad7740dcf
using JLD2

# ╔═╡ 3c07c33b-7f91-4e2d-a27a-09a8924f328c
using DataFrames

# ╔═╡ b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
using CSV

# ╔═╡ 76e2e659-c708-4a08-8263-def6c3dfa930
using Dates

# ╔═╡ 41b0f95e-d806-4167-8199-7392717782f9
using OffsetArrays

# ╔═╡ 342cbcff-5a05-4199-b5b5-6cc0fa70c202
using TemplateMatching

# ╔═╡ f07ab63f-dfb4-493c-a40c-6b3363a858c3
using Plots

# ╔═╡ 356cc4c2-e565-4847-9286-cfc1f835654d
using StatsBase

# ╔═╡ c2f22858-6d3c-4f92-a05b-659deae566f6
md"## Reading Continuous Data"

# ╔═╡ a86b6c8a-9b90-4388-8208-d257ec6951d2
md"""Path of the JLD2 file containing continuous data: $(@bind datapath TextField(default="../data/2021-01-12_20-25-30/data.jld2"))"""

# ╔═╡ f3f17c65-dab6-46da-83a9-87aad1181524
data = load(datapath, "data")

# ╔═╡ f7e0f70d-4069-4039-90b6-d3da559bdcb1
samplefreq = 2 # MHz

# ╔═╡ 57b35604-727b-4eac-bcbd-ea302e4c79a2
md"## Reading events catalogue"

# ╔═╡ b62d7318-3a6d-4442-b6cc-3166b87dff8a
md"""Path of the CSV catalogue $(@bind cataloguepath TextField(default="../data/2021-01-12_20-25-30/cat.csv"))"""

# ╔═╡ d5beeb75-7ca6-4f28-a63b-6e3070fb4c73
begin 
	import Base.+, Base.show
	
	+(x::Dates.UTInstant{Microsecond}, d::Microsecond) = Dates.UTInstant(x.periods + d)

	function show(io::IO, ::MIME"text/plain", x::Dates.UTInstant{Microsecond})
		date = DateTime(Dates.UTM(round(Int, 1e-3 * x.periods.value)))
		show(io, MIME"text/plain"(), date)
	end
end

# ╔═╡ 1bf27aad-b754-4acc-b08d-fe1a86791de8
function DateTimeMicrosecond(y, m, d, h, mi, s, us)
	rata = us + 1_000_000 * (s + 60mi + 3600h + 86400 * Dates.totaldays(y, m, d))
	Dates.UTInstant{Microsecond}(Microsecond(rata))
end

# ╔═╡ 96cf6e3b-1c2b-49f6-979d-0fce28a55a06
starttime = DateTimeMicrosecond(2021, 01, 12, 20, 25, 30, 0)

# ╔═╡ 5d5f3ca7-260d-4dbb-bde1-59f2bd58c0c9
begin
	columns = [:Year, :Month, :Day, :Hour, :Minute, :Second, :North, :East, :Up]
	df = CSV.read(cataloguepath, DataFrame, select=columns)
	sec = floor.(Int, df.Second)
	usec = round.(Int, 1e6 * (df.Second - sec))
	catalogue = DataFrame()
	catalogue.datetime = DateTimeMicrosecond.(df.Year, df.Month, df.Day, df.Hour, df.Minute, sec, usec)
	catalogue.sample = map(x -> round(Int, x.value * samplefreq), catalogue.datetime .- starttime)
	catalogue.north = df.North
	catalogue.east = df.East
	catalogue.up = df.Up
	catalogue
end

# ╔═╡ 6ac6cce3-09ae-421a-b502-b5c8807bf555
md"## Reading Sensors Positions"

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

# ╔═╡ 6ab4a8da-4613-4004-807f-c36876285ebd
md" ## Preparing the templates"

# ╔═╡ 298010b3-c34c-40ba-9e1b-b79f7378f1d4
begin
	displacements_vec = DataFrame[]
	for n in axes(catalogue, 1)
		displacements = DataFrame()
		for s in [:north, :east, :up]
			displacements[!, s] = sensors_positions[!, s] .- catalogue[n, s]
		end
		push!(displacements_vec, displacements)
	end
	displacements_vec
end

# ╔═╡ a1fc5e84-48df-45ea-b6a1-0073ed9f8340
distances_vec = [hypot.(displacements[!, :north], displacements[!, :east], displacements[!, :up]) for displacements in displacements_vec]

# ╔═╡ c7f0582a-8997-4641-b607-0b95a7f73a6c
v_p = 0.67345 # in cm/us

# ╔═╡ 074983c0-2747-4b09-a0c0-4d4adba24c0b
offsets_vec = [@. round(Int, samplefreq * (distances / v_p)) for distances in distances_vec]

# ╔═╡ cd5f912a-fafb-4947-bd70-318151c2c49c
shifts_vec = [catalogue[n, :sample] .+ offsets_vec[n] for n in axes(catalogue, 1)]

# ╔═╡ 033473d3-6958-4e7b-8c95-a9e15d882118
md"""
Template pre:$(@bind t_pre Slider(0:1000, default=100, show_value=true))

Template post:$(@bind t_post Slider(0:1000, default=500, show_value=true))
"""

# ╔═╡ be80dd80-d39d-4338-bb76-c644ea9d970c
cuts_vec = [[shift - t_pre: shift + t_post for shift in shifts] for shifts in shifts_vec]

# ╔═╡ aa23d53f-8333-4c94-8932-fac99ce9196d
templates_vec = [[view(series, cut) for (series, cut) in zip(values(data), cuts)] for cuts in cuts_vec]

# ╔═╡ 368f7fc1-0afc-439f-905a-76af04e8ca8d
md" ## Computing the cross-correlations"

# ╔═╡ 0d122aef-1767-4389-b043-97d07ecef1fe
md"""
tolerance $(@bind tolerance Slider(0:10, default=5, show_value=true))
"""

# ╔═╡ 3cdb79e8-e75a-43cd-a12d-1bffb5a745a0
md"""Path of the JLD2 file containing cross-correlations: $(@bind correlationspath TextField(default="../data/2021-01-12_20-25-30/crosscorrelations.jld2"))"""

# ╔═╡ a7a589fc-4dcc-4137-a423-b8d2e7d41a0e
begin
	if isfile(correlationspath)
		signals_vec = load(correlationspath, "signals_vec")
	else
		data_stream = collect(values(data))
		signals_vec = similar(templates_vec, OffsetVector{Float64, Vector{Float64}})
		Threads.@threads for n = eachindex(signals_vec)
			signals_vec[n] = correlatetemplate(data_stream, templates_vec[n], offsets_vec[n], tolerance)
		end
		save(correlationspath, "signals_vec", signals_vec)
	end
	signals_vec
end

# ╔═╡ abb8b867-39c2-40f2-9d7a-f40bacf8f550
md" ## Finding peaks and Saving Catalogue"

# ╔═╡ 7d1f1903-8342-4715-8e88-3e3d063ea5db
md"""
correlation threshold $(@bind threshold Slider(0.0:0.05:1, default=0.5, show_value=true))

distance (in unit of template length) $(@bind reldistance Slider(0:0.5:5, default=2, show_value=true))
"""

# ╔═╡ 666c8d5f-c1c5-4cb4-9364-bf997b18a683
begin
	peaks_vec = similar(signals_vec, Vector{Int})
	heights_vec = similar(signals_vec, Vector{Float64})
	Threads.@threads for n = eachindex(signals_vec)
		peaks, heights = findpeaks(signals_vec[n], threshold, reldistance * (t_pre + t_post - 1))
		peaks_vec[n] = peaks
		heights_vec[n] = heights
	end
	peaks_vec, heights_vec
end

# ╔═╡ 33bc3279-c65f-47b5-acf4-60853c2bdb91
begin
	templatematch_catalogue = DataFrame()
	for n = eachindex(templates_vec)
		df = DataFrame()
		df.datetime = @. starttime + Microsecond(peaks_vec[n])
		df.sample = peaks_vec[n] .+ t_pre
		df.template .= n
		df.correlation = heights_vec[n]
		df.relative_magnitude = [magnitude(values(data), templates_vec[n], peak .+ offsets_vec[n]) for peak in peaks_vec[n]]
		append!(templatematch_catalogue, df)
	end
	templatematch_catalogue
end

# ╔═╡ 3784a5b2-8621-4bac-87a1-4decc25cf4f5
CSV.write("../data/2021-01-12_20-25-30/templatematch.csv", templatematch_catalogue)

# ╔═╡ 102e4534-6be4-4436-be36-69726a393253
md" ## Statistics of the detections "

# ╔═╡ 45ef95ec-c638-4a94-b631-213cf3170fab
nontemplates = @. !(templatematch_catalogue.correlation > 0.999 && templatematch_catalogue.relative_magnitude ≈ 0.0)

# ╔═╡ feebcd5f-91c4-4a70-a064-e642e3172614
histogram(templatematch_catalogue[nontemplates, :correlation], label=nothing)

# ╔═╡ 476101af-0108-4628-9ac1-f233a456a869
summarystats(templatematch_catalogue[nontemplates, :correlation])

# ╔═╡ 3be237b0-9f19-4e9a-aa13-c6355771b8e8
histogram(templatematch_catalogue[nontemplates, :relative_magnitude], label=nothing)

# ╔═╡ 7453767a-e8f5-4446-b198-fbe48c992aa5
summarystats(templatematch_catalogue[nontemplates, :relative_magnitude])

# ╔═╡ 08f3bb2c-c80e-43ff-814b-e17ae5a912dc
let
	scatter(templatematch_catalogue.sample,
		    catalogue[templatematch_catalogue.template,:sample],
			zcolor=templatematch_catalogue.correlation,
			c=:heat,
			markersize=4exp.(templatematch_catalogue.relative_magnitude),
			ylabel="Template origin [μs]",
			xlabel="Detection origin [μs]",
			colorbar=:right,
			colorbar_title="\nCross-correlation",
   			right_margin = 3Plots.mm,
		    legend=nothing,
			dpi=500)
end

# ╔═╡ Cell order:
# ╠═aa854f02-46bf-4a2b-b1e9-5fb52b99925c
# ╠═8b4426b8-588b-49ad-83c5-a79ed698d704
# ╟─c2f22858-6d3c-4f92-a05b-659deae566f6
# ╟─a86b6c8a-9b90-4388-8208-d257ec6951d2
# ╠═d47c1cf4-4844-4c66-8a21-b62ad7740dcf
# ╠═f3f17c65-dab6-46da-83a9-87aad1181524
# ╠═f7e0f70d-4069-4039-90b6-d3da559bdcb1
# ╟─57b35604-727b-4eac-bcbd-ea302e4c79a2
# ╠═3c07c33b-7f91-4e2d-a27a-09a8924f328c
# ╠═b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
# ╠═76e2e659-c708-4a08-8263-def6c3dfa930
# ╟─b62d7318-3a6d-4442-b6cc-3166b87dff8a
# ╠═1bf27aad-b754-4acc-b08d-fe1a86791de8
# ╠═d5beeb75-7ca6-4f28-a63b-6e3070fb4c73
# ╠═96cf6e3b-1c2b-49f6-979d-0fce28a55a06
# ╠═5d5f3ca7-260d-4dbb-bde1-59f2bd58c0c9
# ╟─6ac6cce3-09ae-421a-b502-b5c8807bf555
# ╠═be090485-7917-464a-8e8c-aec4a2c009d3
# ╟─6ab4a8da-4613-4004-807f-c36876285ebd
# ╠═298010b3-c34c-40ba-9e1b-b79f7378f1d4
# ╠═a1fc5e84-48df-45ea-b6a1-0073ed9f8340
# ╠═c7f0582a-8997-4641-b607-0b95a7f73a6c
# ╠═074983c0-2747-4b09-a0c0-4d4adba24c0b
# ╠═cd5f912a-fafb-4947-bd70-318151c2c49c
# ╠═033473d3-6958-4e7b-8c95-a9e15d882118
# ╠═be80dd80-d39d-4338-bb76-c644ea9d970c
# ╠═aa23d53f-8333-4c94-8932-fac99ce9196d
# ╟─368f7fc1-0afc-439f-905a-76af04e8ca8d
# ╠═41b0f95e-d806-4167-8199-7392717782f9
# ╠═342cbcff-5a05-4199-b5b5-6cc0fa70c202
# ╟─0d122aef-1767-4389-b043-97d07ecef1fe
# ╟─3cdb79e8-e75a-43cd-a12d-1bffb5a745a0
# ╠═a7a589fc-4dcc-4137-a423-b8d2e7d41a0e
# ╟─abb8b867-39c2-40f2-9d7a-f40bacf8f550
# ╟─7d1f1903-8342-4715-8e88-3e3d063ea5db
# ╠═666c8d5f-c1c5-4cb4-9364-bf997b18a683
# ╠═33bc3279-c65f-47b5-acf4-60853c2bdb91
# ╠═3784a5b2-8621-4bac-87a1-4decc25cf4f5
# ╟─102e4534-6be4-4436-be36-69726a393253
# ╠═f07ab63f-dfb4-493c-a40c-6b3363a858c3
# ╠═356cc4c2-e565-4847-9286-cfc1f835654d
# ╠═45ef95ec-c638-4a94-b631-213cf3170fab
# ╟─feebcd5f-91c4-4a70-a064-e642e3172614
# ╟─476101af-0108-4628-9ac1-f233a456a869
# ╟─3be237b0-9f19-4e9a-aa13-c6355771b8e8
# ╟─7453767a-e8f5-4446-b198-fbe48c992aa5
# ╟─08f3bb2c-c80e-43ff-814b-e17ae5a912dc
