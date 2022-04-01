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

# ╔═╡ 7ddcab3c-a100-44be-baa3-af28cc08b2d7
using Interpolations

# ╔═╡ 341b9066-dcb3-4d99-b985-82e24c11f28e
using Optim

# ╔═╡ c2f22858-6d3c-4f92-a05b-659deae566f6
md"## Reading Continuous Data"

# ╔═╡ a86b6c8a-9b90-4388-8208-d257ec6951d2
md"""Path of the JLD2 file containing continuous data: $(@bind datapath TextField(default="../data/2021-01-12_20-25-30/data.jld2"))"""

# ╔═╡ f3f17c65-dab6-46da-83a9-87aad1181524
data = load(datapath, "data")

# ╔═╡ 8271f7e6-d8ac-472a-9538-75d0afacf42e
samplefreq = 2 # MHz

# ╔═╡ 57b35604-727b-4eac-bcbd-ea302e4c79a2
md"## Reading events catalogue"

# ╔═╡ b62d7318-3a6d-4442-b6cc-3166b87dff8a
md"""Path of the CSV catalogue $(@bind cataloguepath TextField(default="../data/2021-01-12_20-25-30/cat.csv"))"""

# ╔═╡ f3139f19-e386-4bb7-8aca-fa54062d4093
begin 
	import Base.+, Base.show
	
	+(x::Dates.UTInstant{Microsecond}, d::Microsecond) = Dates.UTInstant(x.periods + d)

	function show(io::IO, ::MIME"text/plain", x::Dates.UTInstant{Microsecond})
		date = DateTime(Dates.UTM(round(Int, 1e-3 * x.periods.value)))
		show(io, MIME"text/plain"(), date)
	end
end

# ╔═╡ 9de113e1-36ce-4251-9fa5-5baa27941a85
function DateTimeMicrosecond(y, m, d, h, mi, s, us)
	rata = us + 1_000_000 * (s + 60mi + 3600h + 86400 * Dates.totaldays(y, m, d))
	Dates.UTInstant{Microsecond}(Microsecond(rata))
end

# ╔═╡ 96cf6e3b-1c2b-49f6-979d-0fce28a55a06
starttime = DateTimeMicrosecond(2021, 01, 12, 20, 25, 30, 0)

# ╔═╡ 5d5f3ca7-260d-4dbb-bde1-59f2bd58c0c9
begin
	df = CSV.read(cataloguepath, DataFrame)
	sec = floor.(Int, df.Second)
	usec = round.(Int, 1e6 * (df.Second - sec))
	templates = DataFrame()
	templates.datetime = DateTimeMicrosecond.(df.Year, df.Month, df.Day, df.Hour, df.Minute, sec, usec)
	templates.sample = map(x -> round(Int, x.value * samplefreq), templates.datetime .- starttime)
	templates.north = df.North
	templates.east = df.East
	templates.up = df.Up
	templates
end

# ╔═╡ 6ac6cce3-09ae-421a-b502-b5c8807bf555
md"## Reading Sensors Positions"

# ╔═╡ be090485-7917-464a-8e8c-aec4a2c009d3
begin
	sensors_positions = CSV.read("../data/2021-01-12_20-25-30/passive-xyz.csv", DataFrame, header=[:north, :east, :up])
	sensors_positions .*= 100
	sensors_positions.sensor = axes(sensors_positions, 1) .- 1
	filter!(r -> r.sensor in keys(data), sensors_positions)
	sensors_positions
end

# ╔═╡ 6ab4a8da-4613-4004-807f-c36876285ebd
md" ## Preparing a single template"

# ╔═╡ fe42e7e6-2c56-4f11-aae6-ec06ad7a219f
md"""template number: $(@bind template_num Select(1:nrow(templates), default=116))"""

# ╔═╡ 298010b3-c34c-40ba-9e1b-b79f7378f1d4
begin
	displacement = DataFrame()
	for s in [:north, :east, :up]
		displacement[!, s] = sensors_positions[!, s] .- templates[template_num, s]
	end
	displacement
end

# ╔═╡ a1fc5e84-48df-45ea-b6a1-0073ed9f8340
distances = hypot.(displacement[!, :north], displacement[!, :east], displacement[!, :up]) # in cm

# ╔═╡ c7f0582a-8997-4641-b607-0b95a7f73a6c
v_p = 0.67345 # in cm/us

# ╔═╡ 074983c0-2747-4b09-a0c0-4d4adba24c0b
offsets = @. round(Int, samplefreq * (distances / v_p))

# ╔═╡ cd5f912a-fafb-4947-bd70-318151c2c49c
shifts = templates[template_num, :sample] .+ offsets

# ╔═╡ 033473d3-6958-4e7b-8c95-a9e15d882118
md"""
Template window pre (in samples): $(@bind t_pre Slider(0:10:1000, default=100, show_value=true))

Template window post (in samples): $(@bind t_post Slider(0:10:1000, default=500, show_value=true))
"""

# ╔═╡ be80dd80-d39d-4338-bb76-c644ea9d970c
cuts = [shift - t_pre: shift + t_post for shift in shifts]

# ╔═╡ aa23d53f-8333-4c94-8932-fac99ce9196d
template_data = [view(series, cut) for (series, cut) in zip(values(data), cuts)]

# ╔═╡ 4ef908f4-d8dc-4159-97ad-fc8d4621a652
let
	template_aligned = OffsetVector.(template_data, shifts)
	xlim = minimum(first(axes(series, 1)) for series in template_aligned), maximum(last(axes(series, 1)) for series in template_aligned)
	plots = []
	for (nch, series) = zip(keys(data), template_aligned)
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
	numchannels = length(template_data)
	plot(plots..., layout=(numchannels, 1), size=(800, 100numchannels), dpi=300)
end

# ╔═╡ 368f7fc1-0afc-439f-905a-76af04e8ca8d
md" ## Computing the cross-correlation"

# ╔═╡ 0d122aef-1767-4389-b043-97d07ecef1fe
md"""
tolerance $(@bind tolerance Slider(0:10, default=5, show_value=true))
"""

# ╔═╡ a7a589fc-4dcc-4137-a423-b8d2e7d41a0e
signal = correlatetemplate(collect(values(data)), template_data, offsets, tolerance)

# ╔═╡ 91924a12-6ee4-4d1b-8c49-9b54142cf733
let
	template_ref = templates[template_num, :sample] - t_pre
	window = t_pre + t_post - 1
	a = max(first(axes(signal, 1)), template_ref - div(window, 2))
	b = min(last(axes(signal, 1)), template_ref + div(window, 2))
	plot(a:b, signal[a:b], label=nothing, dpi=300)
end

# ╔═╡ abb8b867-39c2-40f2-9d7a-f40bacf8f550
md" ## Finding peaks"

# ╔═╡ 7d1f1903-8342-4715-8e88-3e3d063ea5db
md"""
correlation threshold $(@bind threshold Slider(0.0:0.05:1, default=0.4, show_value=true))

distance (in unit of template length) $(@bind reldistance Slider(0:0.5:5, default=2, show_value=true))
"""

# ╔═╡ 666c8d5f-c1c5-4cb4-9364-bf997b18a683
peaks, heights = findpeaks(signal, threshold, reldistance * (t_pre + t_post - 1))

# ╔═╡ 96deaa3d-fc40-466f-bf1b-d8d3dd58dcd4
md"""
Peak number $(@bind peak_num Slider(axes(peaks, 1), default=1, show_value=true)) / $(length(peaks))
"""

# ╔═╡ cc0c3e0a-4208-4339-bea1-75396b2778bd
let 
	plt = plot(ylims=(0.0, 1.0))
	scale = round(Int, reldistance * (t_post + t_pre - 1))
	a = max(first(axes(signal, 1)), peaks[peak_num] - scale)
	b = min(last(axes(signal, 1)), peaks[peak_num] + scale)
	plot!(plt, a:b, signal[a:b], label=nothing)
	plot!(plt, peaks[a .< peaks .< b], heights[a.< peaks .< b], 
		seriestype=:scatter, label=nothing, dpi=300)
	plt
end

# ╔═╡ b6848d59-e97e-4662-a2ee-78612694d6d8
md"## Relocalization"

# ╔═╡ a899193d-8685-48b3-a3ff-d91973bd6834
function find_toa(trace, waveform, center, left, tol)
	right = size(waveform, 1) - left - 1
	a, b = center - left - tol, center + right + tol
	cc = OffsetVector(crosscorrelate(view(trace, a:b), waveform), -tol - 1)
	cc_max, n_max = findmax(cc)
	if firstindex(cc) < n_max < lastindex(cc) && cc_max < 1.0 && cc_max ≉ 1.0
		y1, y2, y3 = cc[n_max - 1:n_max + 1]
		delta = (y1 - y3) / 2(y1 - 2y2 + y3)
		center + n_max + delta, min(1.0, y2 + (y3 - y1) / 4delta)
	else
		center + Float64(n_max), cc_max
	end
end

# ╔═╡ bd5e4f6f-c7e3-440c-9f02-fb3613625508
toas_corrs = find_toa.(collect(values(data)), template_data, (peaks[peak_num] + t_pre) .+ offsets, t_pre, tolerance)

# ╔═╡ dd0e7c9c-15d7-4e50-aff2-dca482d25220
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

# ╔═╡ 6d6d27c9-3b1a-4646-8420-61a775f7842f
let 
	peak = peaks[peak_num]
	peak_start = peak - div(t_pre, 2)
	peak_stop = peak + (t_pre + t_post - 1) + div(t_post, 2)
	data_start = minimum(series -> first(axes(series, 1)), values(data))
	data_stop = maximum(series -> last(axes(series, 1)), values(data))
	xlim = max(data_start, peak_start), min(data_stop, peak_stop)
	
	plots = []
	for (n, channel_num) = enumerate(keys(data))
		series_start = 	max(first(axes(data[channel_num], 1)), peak_start)
		series_stop = min(last(axes(data[channel_num], 1)), peak_stop)
		data_series = OffsetVector(view(data[channel_num], series_start:series_stop), series_start:series_stop)
		
		toa_f, toa_i = modf(toas_corrs[n][1] - t_pre)
		template_interpolated = subsampleshift(template_data[n], toa_f)
		template_aligned = OffsetVector(template_interpolated, Int(toa_i))
		
		data_inf, data_sup = extrema(data_series)
		data_center = mean(data_series)
		template_inf, template_sup = extrema(template_interpolated)
		template_center = mean(template_aligned)
		β = (data_sup - data_inf) / (template_sup - template_inf)
		push!(plots,
			  plot([axes(data_series, 1), axes(template_aligned, 1)], 
				   [data_series, data_center .+ β .* (template_aligned .- template_center)];
				   xlim,
				   ylim=(data_inf, data_sup),
				   title="Channel $channel_num",
				   titlealign=:left,
				   showaxis=false,
				   ticks=false,
				   label=nothing))
	end
	numchannels = length(template_data)
	plot(plots..., layout=(numchannels, 1), size=(800, 100numchannels), dpi=300)
end

# ╔═╡ 2911437d-4775-4a1b-819e-adf6e91acd61
measurements = let
	measurements = copy(sensors_positions)
	measurements.toas = [sample / samplefreq for (sample, _) in toas_corrs]
	measurements.cc = [cc for (_, cc) in toas_corrs]
	filter(row -> row.cc > 0.5, measurements)
end

# ╔═╡ abbf04cf-ccaf-45a1-8f3e-95aa770d6f4b
function residue(s)
	x, t = s[1:3], s[4]
	Δx_squared = (transpose(Vector(x)) .- measurements[!, [:north, :east, :up]]) .^ 2
	d_squared =  (v_p .* (t .- measurements.toas)) .^ 2
	sum((sum(eachcol(Δx_squared)) .- d_squared) .^ 2)
end

# ╔═╡ 2e6b7b35-9512-406d-9a97-f1f94bea8732
if nrow(measurements) > 4
	x0 = Vector(templates[template_num, [:north, :east, :up]])
	t0 = (peaks[peak_num] + t_pre) / samplefreq
	candidate = optimize(residue, [x0; t0])
end

# ╔═╡ 32ab23a5-60b6-4f07-8c0f-cfa2e0a9d7e1
candidate.minimizer

# ╔═╡ 102e4534-6be4-4436-be36-69726a393253
md" ## Statistics of the detections "

# ╔═╡ 33bc3279-c65f-47b5-acf4-60853c2bdb91
begin
	catalogue = DataFrame()
	samples = peaks .+ t_pre
	catalogue.datetimes = @. starttime + Microsecond(samples)
	catalogue.sample = samples
	catalogue.correlation = heights
	catalogue.rel_mag = [magnitude(values(data), template_data, peak .+ offsets) for peak in peaks]
	catalogue
end

# ╔═╡ 471f21ba-23d3-42f4-8d62-12fe2f36124f
histogram(collect(signal), label=nothing)

# ╔═╡ 6e2e9172-4cf5-45f3-a239-12c09b6ca3da
summarystats(collect(signal))

# ╔═╡ feebcd5f-91c4-4a70-a064-e642e3172614
histogram(catalogue[!, :correlation], label=nothing)

# ╔═╡ 476101af-0108-4628-9ac1-f233a456a869
summarystats(catalogue[!, :correlation])

# ╔═╡ 3be237b0-9f19-4e9a-aa13-c6355771b8e8
histogram(catalogue[!, :rel_mag], label=nothing)

# ╔═╡ Cell order:
# ╠═aa854f02-46bf-4a2b-b1e9-5fb52b99925c
# ╟─c2f22858-6d3c-4f92-a05b-659deae566f6
# ╠═8b4426b8-588b-49ad-83c5-a79ed698d704
# ╠═d47c1cf4-4844-4c66-8a21-b62ad7740dcf
# ╟─a86b6c8a-9b90-4388-8208-d257ec6951d2
# ╠═f3f17c65-dab6-46da-83a9-87aad1181524
# ╠═8271f7e6-d8ac-472a-9538-75d0afacf42e
# ╟─57b35604-727b-4eac-bcbd-ea302e4c79a2
# ╠═3c07c33b-7f91-4e2d-a27a-09a8924f328c
# ╠═b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
# ╠═76e2e659-c708-4a08-8263-def6c3dfa930
# ╟─b62d7318-3a6d-4442-b6cc-3166b87dff8a
# ╠═9de113e1-36ce-4251-9fa5-5baa27941a85
# ╠═96cf6e3b-1c2b-49f6-979d-0fce28a55a06
# ╠═f3139f19-e386-4bb7-8aca-fa54062d4093
# ╠═5d5f3ca7-260d-4dbb-bde1-59f2bd58c0c9
# ╟─6ac6cce3-09ae-421a-b502-b5c8807bf555
# ╠═be090485-7917-464a-8e8c-aec4a2c009d3
# ╟─6ab4a8da-4613-4004-807f-c36876285ebd
# ╠═41b0f95e-d806-4167-8199-7392717782f9
# ╠═95f2bd0b-7d8b-43b9-b2f9-cfa6628ec485
# ╠═356cc4c2-e565-4847-9286-cfc1f835654d
# ╠═730340bb-8176-4282-8559-8d60b1d28d57
# ╠═f07ab63f-dfb4-493c-a40c-6b3363a858c3
# ╟─fe42e7e6-2c56-4f11-aae6-ec06ad7a219f
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
# ╟─91924a12-6ee4-4d1b-8c49-9b54142cf733
# ╟─abb8b867-39c2-40f2-9d7a-f40bacf8f550
# ╟─7d1f1903-8342-4715-8e88-3e3d063ea5db
# ╠═666c8d5f-c1c5-4cb4-9364-bf997b18a683
# ╟─96deaa3d-fc40-466f-bf1b-d8d3dd58dcd4
# ╟─cc0c3e0a-4208-4339-bea1-75396b2778bd
# ╟─b6848d59-e97e-4662-a2ee-78612694d6d8
# ╠═7ddcab3c-a100-44be-baa3-af28cc08b2d7
# ╠═341b9066-dcb3-4d99-b985-82e24c11f28e
# ╠═a899193d-8685-48b3-a3ff-d91973bd6834
# ╠═bd5e4f6f-c7e3-440c-9f02-fb3613625508
# ╠═dd0e7c9c-15d7-4e50-aff2-dca482d25220
# ╟─6d6d27c9-3b1a-4646-8420-61a775f7842f
# ╠═2911437d-4775-4a1b-819e-adf6e91acd61
# ╠═abbf04cf-ccaf-45a1-8f3e-95aa770d6f4b
# ╠═2e6b7b35-9512-406d-9a97-f1f94bea8732
# ╠═32ab23a5-60b6-4f07-8c0f-cfa2e0a9d7e1
# ╟─102e4534-6be4-4436-be36-69726a393253
# ╠═33bc3279-c65f-47b5-acf4-60853c2bdb91
# ╠═471f21ba-23d3-42f4-8d62-12fe2f36124f
# ╠═6e2e9172-4cf5-45f3-a239-12c09b6ca3da
# ╠═feebcd5f-91c4-4a70-a064-e642e3172614
# ╠═476101af-0108-4628-9ac1-f233a456a869
# ╠═3be237b0-9f19-4e9a-aa13-c6355771b8e8
