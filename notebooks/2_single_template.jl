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

# ╔═╡ 95f2bd0b-7d8b-43b9-b2f9-cfa6628ec485
using LinearAlgebra

# ╔═╡ 356cc4c2-e565-4847-9286-cfc1f835654d
using StatsBase

# ╔═╡ 730340bb-8176-4282-8559-8d60b1d28d57
using Statistics

# ╔═╡ f07ab63f-dfb4-493c-a40c-6b3363a858c3
using Plots

# ╔═╡ 342cbcff-5a05-4199-b5b5-6cc0fa70c202
using TemplateMatching

# ╔═╡ c2f22858-6d3c-4f92-a05b-659deae566f6
md"## Reading Continuous Data"

# ╔═╡ a86b6c8a-9b90-4388-8208-d257ec6951d2
md"""Path of the JLD2 file containing continuous data: $(@bind datapath TextField(default="../data/2021-01-12_20-25-30/data.jld2"))"""

# ╔═╡ f3f17c65-dab6-46da-83a9-87aad1181524
data = load(datapath, "data")

# ╔═╡ 57b35604-727b-4eac-bcbd-ea302e4c79a2
md"## Reading events catalogue"

# ╔═╡ b62d7318-3a6d-4442-b6cc-3166b87dff8a
md"""Path of the CSV catalogue $(@bind cataloguepath TextField(default="../data/2021-01-12_20-25-30/cat2021-01-12_20-25-30_gab6_.2021.5.10.23.37.11.3267.csv"))"""

# ╔═╡ 96cf6e3b-1c2b-49f6-979d-0fce28a55a06
starttime = DateTime(2021, 01, 12, 20, 25, 30)

# ╔═╡ 5d5f3ca7-260d-4dbb-bde1-59f2bd58c0c9
begin
	columns = [:Year, :Month, :Day, :Hour, :Minute, :Second, :North, :East, :Up]
	catalogue = CSV.read(cataloguepath, DataFrame, select=columns)
	sec = floor.(Int, catalogue.Second)
	usec = round.(Int, 1e6 * (catalogue.Second - sec))
	instants = @. usec + 1_000_000 * (sec + 60 * catalogue.Minute + 3600 * catalogue.Hour + 86_400 * Dates.totaldays(catalogue.Year, catalogue.Month, catalogue.Day))
	catalogue.sample = @. instants - Microsecond(starttime.instant.periods).value
	catalogue
end

# ╔═╡ a816019d-529e-44cd-b0da-e6985fadf494
catalogue[!, [:sample, :North, :East, :Up]]

# ╔═╡ 6ac6cce3-09ae-421a-b502-b5c8807bf555
md"## Reading Sensors Positions"

# ╔═╡ be090485-7917-464a-8e8c-aec4a2c009d3
begin
	sensors_positions = CSV.read("../data/2021-01-12_20-25-30/passive-xyz.csv", DataFrame, header=[:North, :East, :Up])
	for s = [:North, :East, :Up]
		sensors_positions[!, s] .*= 100
	end
	sensors_positions.sensor = axes(sensors_positions, 1) .- 1
	filter!(r -> r.sensor in keys(data), sensors_positions)
	sensors_positions
end

# ╔═╡ 6ab4a8da-4613-4004-807f-c36876285ebd
md" ## Preparing a single template"

# ╔═╡ fe42e7e6-2c56-4f11-aae6-ec06ad7a219f
md"""template number: $(@bind template_num Select(1:nrow(catalogue), default=116))"""

# ╔═╡ 726a60e7-9b74-4c6e-afdc-abdc151d9583
event_coordinates = Vector(catalogue[template_num, [:North, :East, :Up]])

# ╔═╡ 298010b3-c34c-40ba-9e1b-b79f7378f1d4
begin
	displacement = DataFrame()
	for s in [:North, :East, :Up]
		displacement[!, s] = sensors_positions[!, s] .- catalogue[template_num, s]
	end
	displacement
end

# ╔═╡ a1fc5e84-48df-45ea-b6a1-0073ed9f8340
distances = hypot.(displacement[!, :North], displacement[!, :East], displacement[!, :Up])

# ╔═╡ c7f0582a-8997-4641-b607-0b95a7f73a6c
v_p = 0.67345 # in cm/us

# ╔═╡ 074983c0-2747-4b09-a0c0-4d4adba24c0b
offsets = @. round(Int, distances / v_p)

# ╔═╡ cd5f912a-fafb-4947-bd70-318151c2c49c
shifts = catalogue[template_num, :sample] .+ offsets

# ╔═╡ 033473d3-6958-4e7b-8c95-a9e15d882118
md"""
Template pre:$(@bind t_pre Slider(0:400, default=200, show_value=true))

Template post:$(@bind t_post Slider(0:400, default=200, show_value=true))
"""

# ╔═╡ be80dd80-d39d-4338-bb76-c644ea9d970c
cuts = [shift - t_pre: shift + t_post for shift in shifts]

# ╔═╡ aa23d53f-8333-4c94-8932-fac99ce9196d
template = [view(series, cut) for (series, cut) in zip(values(data), cuts)]

# ╔═╡ 4ef908f4-d8dc-4159-97ad-fc8d4621a652
let
	xlim = minimum(first(axes(series, 1)) for series in template), maximum(last(axes(series, 1)) for series in template)
	plots = []
	for (nch, series, shift, cut) = zip(keys(data), template, shifts, cuts)
		ylim = minimum(series), maximum(series)
		push!(plots, 
			  plot(axes(series, 1), series;
				   xlim,
				   ylim,
			       title="Channel $nch",
				   titlealign=:left,
				   showaxis=false,
				   ticks=false,
				   label=nothing))
	end
	numchannels = length(template)
	plot(plots..., layout=(numchannels, 1), size=(800, 100numchannels), dpi=300)
end

# ╔═╡ 368f7fc1-0afc-439f-905a-76af04e8ca8d
md" ## Computing the cross-correlation"

# ╔═╡ 0d122aef-1767-4389-b043-97d07ecef1fe
md"""
tolerance $(@bind tolerance Slider(0:10, default=5, show_value=true))
"""

# ╔═╡ a7a589fc-4dcc-4137-a423-b8d2e7d41a0e
signal = correlatetemplate(collect(values(data)), template, offsets, tolerance)

# ╔═╡ 91924a12-6ee4-4d1b-8c49-9b54142cf733
let
	a = max(first(axes(signal, 1)), catalogue[template_num, :sample] - 2t_pre)
	b = min(last(axes(signal, 1)), catalogue[template_num, :sample] + t_post - t_pre)
	plot(a:b, signal[a:b], label=nothing, dpi=300)
end

# ╔═╡ abb8b867-39c2-40f2-9d7a-f40bacf8f550
md" ## Finding peaks"

# ╔═╡ 7d1f1903-8342-4715-8e88-3e3d063ea5db
md"""
correlation threshold $(@bind threshold Slider(0.0:0.05:1, default=0.5, show_value=true))

distance (in unit of template length) $(@bind reldistance Slider(0:0.5:5, default=2, show_value=true))
"""

# ╔═╡ 666c8d5f-c1c5-4cb4-9364-bf997b18a683
peaks, heights = findpeaks(signal, threshold, reldistance * (t_pre + t_post - 1))

# ╔═╡ 96deaa3d-fc40-466f-bf1b-d8d3dd58dcd4
md"""
Peak number $(@bind window_pos Slider(axes(peaks, 1), default=1, show_value=true))
"""

# ╔═╡ cc0c3e0a-4208-4339-bea1-75396b2778bd
let plt = plot(ylims=(0.0, 1.0)), scale = round(Int, reldistance * (t_post + t_pre))
	a = max(first(axes(signal, 1)), peaks[window_pos] - scale)
	b = min(last(axes(signal, 1)), peaks[window_pos] + scale)
	plot!(plt, a:b, signal[a:b], label=nothing)
	plot!(plt, peaks[a .< peaks .< b], heights[a.< peaks .< b], 
		seriestype=:scatter, label=nothing, dpi=300)
	plt
end

# ╔═╡ ff450e8c-83ee-4e70-b9f4-31d6e1e99b62
catalogue[template_num, [:Year, :Month, :Day, :Hour, :Minute, :Second, :sample]]

# ╔═╡ 33bc3279-c65f-47b5-acf4-60853c2bdb91
begin
	datetimes = @. starttime + Microsecond(peaks - t_pre)
	output_catalogue = DataFrame()
	for (r, s) in [(:year, :Year), (:month, :Month), (:day, :Day), (:hour, :Hour), (:minute, :Minute)]
		output_catalogue[!, s] = getproperty(Dates, r).(datetimes)
	end
	output_catalogue.Second = @. Dates.second(datetimes) + 1e-3 * Dates.millisecond(datetimes)
	output_catalogue.sample = peaks .- t_pre
	output_catalogue.template .= template_num
	output_catalogue.correlation = heights
	output_catalogue.relative_magnitude = [magnitude(values(data), template, peak .+ offsets .- t_pre) for peak in peaks]
	output_catalogue
end

# ╔═╡ 102e4534-6be4-4436-be36-69726a393253
md" ## Statistics of the detections "

# ╔═╡ 471f21ba-23d3-42f4-8d62-12fe2f36124f
histogram(collect(signal), label=nothing)

# ╔═╡ 6e2e9172-4cf5-45f3-a239-12c09b6ca3da
summarystats(collect(signal))

# ╔═╡ feebcd5f-91c4-4a70-a064-e642e3172614
histogram(output_catalogue[!, :correlation], label=nothing)

# ╔═╡ 476101af-0108-4628-9ac1-f233a456a869
summarystats(output_catalogue[!, :correlation])

# ╔═╡ 3be237b0-9f19-4e9a-aa13-c6355771b8e8
histogram(output_catalogue[!, :relative_magnitude], label=nothing)

# ╔═╡ Cell order:
# ╠═aa854f02-46bf-4a2b-b1e9-5fb52b99925c
# ╟─c2f22858-6d3c-4f92-a05b-659deae566f6
# ╠═8b4426b8-588b-49ad-83c5-a79ed698d704
# ╟─a86b6c8a-9b90-4388-8208-d257ec6951d2
# ╠═d47c1cf4-4844-4c66-8a21-b62ad7740dcf
# ╠═f3f17c65-dab6-46da-83a9-87aad1181524
# ╟─57b35604-727b-4eac-bcbd-ea302e4c79a2
# ╠═3c07c33b-7f91-4e2d-a27a-09a8924f328c
# ╠═b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
# ╠═76e2e659-c708-4a08-8263-def6c3dfa930
# ╟─b62d7318-3a6d-4442-b6cc-3166b87dff8a
# ╠═96cf6e3b-1c2b-49f6-979d-0fce28a55a06
# ╠═5d5f3ca7-260d-4dbb-bde1-59f2bd58c0c9
# ╠═a816019d-529e-44cd-b0da-e6985fadf494
# ╟─6ac6cce3-09ae-421a-b502-b5c8807bf555
# ╠═be090485-7917-464a-8e8c-aec4a2c009d3
# ╟─6ab4a8da-4613-4004-807f-c36876285ebd
# ╠═41b0f95e-d806-4167-8199-7392717782f9
# ╠═95f2bd0b-7d8b-43b9-b2f9-cfa6628ec485
# ╠═356cc4c2-e565-4847-9286-cfc1f835654d
# ╟─fe42e7e6-2c56-4f11-aae6-ec06ad7a219f
# ╠═726a60e7-9b74-4c6e-afdc-abdc151d9583
# ╠═730340bb-8176-4282-8559-8d60b1d28d57
# ╠═f07ab63f-dfb4-493c-a40c-6b3363a858c3
# ╠═298010b3-c34c-40ba-9e1b-b79f7378f1d4
# ╠═a1fc5e84-48df-45ea-b6a1-0073ed9f8340
# ╠═c7f0582a-8997-4641-b607-0b95a7f73a6c
# ╠═074983c0-2747-4b09-a0c0-4d4adba24c0b
# ╠═cd5f912a-fafb-4947-bd70-318151c2c49c
# ╟─033473d3-6958-4e7b-8c95-a9e15d882118
# ╠═be80dd80-d39d-4338-bb76-c644ea9d970c
# ╠═aa23d53f-8333-4c94-8932-fac99ce9196d
# ╟─4ef908f4-d8dc-4159-97ad-fc8d4621a652
# ╟─368f7fc1-0afc-439f-905a-76af04e8ca8d
# ╠═342cbcff-5a05-4199-b5b5-6cc0fa70c202
# ╟─0d122aef-1767-4389-b043-97d07ecef1fe
# ╠═a7a589fc-4dcc-4137-a423-b8d2e7d41a0e
# ╠═91924a12-6ee4-4d1b-8c49-9b54142cf733
# ╟─abb8b867-39c2-40f2-9d7a-f40bacf8f550
# ╟─7d1f1903-8342-4715-8e88-3e3d063ea5db
# ╠═666c8d5f-c1c5-4cb4-9364-bf997b18a683
# ╟─96deaa3d-fc40-466f-bf1b-d8d3dd58dcd4
# ╠═cc0c3e0a-4208-4339-bea1-75396b2778bd
# ╠═ff450e8c-83ee-4e70-b9f4-31d6e1e99b62
# ╠═33bc3279-c65f-47b5-acf4-60853c2bdb91
# ╟─102e4534-6be4-4436-be36-69726a393253
# ╠═471f21ba-23d3-42f4-8d62-12fe2f36124f
# ╠═6e2e9172-4cf5-45f3-a239-12c09b6ca3da
# ╠═feebcd5f-91c4-4a70-a064-e642e3172614
# ╠═476101af-0108-4628-9ac1-f233a456a869
# ╠═3be237b0-9f19-4e9a-aa13-c6355771b8e8
