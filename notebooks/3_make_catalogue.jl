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

# ╔═╡ 4249c937-1b33-4aca-8351-d98a19c584fa
using Tables

# ╔═╡ 87df9037-d497-44e8-836d-3c34a58c5d07
using StatsBase

# ╔═╡ c2f22858-6d3c-4f92-a05b-659deae566f6
md"## Reading Continuous Data"

# ╔═╡ a86b6c8a-9b90-4388-8208-d257ec6951d2
md"""Path of the JLD2 file containing continuous data: $(@bind datapath TextField(default="../data/2021-01-12_20-25-30/data.jld2"))"""

# ╔═╡ f3f17c65-dab6-46da-83a9-87aad1181524
data = load(datapath, "data")

# ╔═╡ f7e0f70d-4069-4039-90b6-d3da559bdcb1
samplefreq = 2 # MHz

# ╔═╡ c7f0582a-8997-4641-b607-0b95a7f73a6c
v_p = 0.67345 # in cm/us

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
catalogue = let
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

# ╔═╡ cfceae8c-0967-4a47-93f5-2b6e3073624c
md"""Path of the CSV of sensors positions $(@bind sensorspospath TextField(default="../data/2021-01-12_20-25-30/passive-xyz.csv"))"""

# ╔═╡ be090485-7917-464a-8e8c-aec4a2c009d3
sensors_positions = let
	sensors_positions = CSV.read(sensorspospath, DataFrame, header=[:north, :east, :up])
	for s = [:north, :east, :up]
		sensors_positions[!, s] .*= 100
	end
	sensors_positions.sensor = axes(sensors_positions, 1) .- 1
	filter!(r -> r.sensor in keys(data), sensors_positions)
	sensors_positions
end

# ╔═╡ 368f7fc1-0afc-439f-905a-76af04e8ca8d
md"## Processing the catalogue"

# ╔═╡ 033473d3-6958-4e7b-8c95-a9e15d882118
md"""
Template pre:$(@bind t_pre Slider(0:1000, default=100, show_value=true))
"""

# ╔═╡ 1ba147d6-2609-4d83-a734-6bf91445a23c
md"""Template post:$(@bind t_post Slider(0:1000, default=500, show_value=true))
"""

# ╔═╡ 7d1f1903-8342-4715-8e88-3e3d063ea5db
md"""
Signal height threshold $(@bind height_threshold Slider(0.0:0.05:1, default=0.3, show_value=true))
"""

# ╔═╡ 8790874d-322e-4aea-a60b-4717e653afc9
md"""Distance (in unit of template length) $(@bind reldistance Slider(0:0.5:5, default=2, show_value=true))"""

# ╔═╡ 0d122aef-1767-4389-b043-97d07ecef1fe
md"""
Tolerance $(@bind tolerance Slider(0:10, default=5, show_value=true))
"""

# ╔═╡ ef346f24-ab46-407b-95a0-2bce10b41ac8
md"""Cross-correlation threshold $(@bind cc_threshold Slider(0.0:0.05:1, default=0.5, show_value=true))"""

# ╔═╡ c4b39f91-8f0d-4f54-90b4-3c035a5b125f
md"""Minimum number of channels above threshold $(@bind nch_min Slider(0:length(data), default=4, show_value=true))"""

# ╔═╡ afe1c12f-23d3-47af-bfe8-b59ca7018593
function process_match(data_stream, template, match_offsets, guess, correlation_threshold, nch_min)
	toas = estimatetoa.(data_stream, template, match_offsets, tolerance)
	readings = copy(sensors_positions)
	readings.index = 1:nrow(readings)
	readings.toas = [(sample + t_pre) / samplefreq for (sample, _) in toas]
	readings.cc = [cc for (_, cc) in toas]
	filter!(row -> row.cc > correlation_threshold, readings)
	crosscorrelation = mean(readings.cc)
	relative_magnitude = magnitude(data_stream[readings.index], template[readings.index], match_offsets[readings.index])
	nch = nrow(readings)
	if nch >= nch_min
		readings_matrix = eachrow(Matrix(readings[!, [:north, :east, :up, :toas]]))
		candidate = locate(readings_matrix, v_p, guess)
		north, east, up, origin_time = candidate.minimizer
		multilateration_residual = candidate.minimum
	else
		north, east, up, origin_time = guess
		multilateration_residual = missing
	end
	(; north, east, up, origin_time, relative_magnitude, crosscorrelation, nch, multilateration_residual)
end

# ╔═╡ 666c8d5f-c1c5-4cb4-9364-bf997b18a683
augmented_catalogue = let
	data_stream = collect(values(data))
	matches_vec = Vector{DataFrame}(undef, nrow(catalogue))
	distance_threshold = reldistance * (t_pre + t_post - 1)
	Threads.@threads for n = eachindex(matches_vec)
		# Preparing template
		displacement = DataFrame()
		for s in [:north, :east, :up]
			displacement[!, s] = sensors_positions[!, s] .- catalogue[n, s]
		end
		distances = map(Base.splat(hypot), eachrow(displacement[!, [:north, :east, :up]]))
		offsets = @. round(Int, samplefreq * (distances / v_p))
		shifts = catalogue[n, :sample] .+ offsets
		cuts = [shift - t_pre: shift + t_post for shift in shifts]
		template_data = [view(series, cut) for (series, cut) in zip(data_stream, cuts)]
		# Computing cross-correlation and finding peaks
		signal = correlatetemplate(data_stream, template_data, offsets, tolerance)
		peaks, heights = findpeaks(signal, height_threshold, distance_threshold)
		# Computing magnitudes and relocating events
		matches = DataFrame()
		matches.peak_sample = peaks .+ t_pre
		matches.peak_height = heights
		matches.template .= n
		x0 = Vector(catalogue[n, [:north, :east, :up]])
		matches_data = [process_match(data_stream, 
			                          template_data, 
			                          peak .+ offsets, 
			                          [x0; (peak + t_pre) / samplefreq],
			                          cc_threshold, 
			                          nch_min) 
			            for peak in peaks]
		matches = hcat(matches, DataFrame(matches_data))
		matches_vec[n] = matches
	end
	augmented_catalogue = reduce(vcat, matches_vec)
	tokeep = TemplateMatching.selectbypeaksdistance(augmented_catalogue.origin_time, augmented_catalogue.crosscorrelation .* augmented_catalogue.nch, distance_threshold / samplefreq)
	augmented_catalogue[tokeep, :]
end

# ╔═╡ b50df30e-f7c3-44f2-b759-3e4ca31003d9
md"""Path where to output the augmented catalogue $(@bind outputpath TextField(default="../data/2021-01-12_20-25-30/augmented_catalogue.csv"))"""

# ╔═╡ 3784a5b2-8621-4bac-87a1-4decc25cf4f5
CSV.write(outputpath, augmented_catalogue)

# ╔═╡ Cell order:
# ╠═aa854f02-46bf-4a2b-b1e9-5fb52b99925c
# ╠═8b4426b8-588b-49ad-83c5-a79ed698d704
# ╟─c2f22858-6d3c-4f92-a05b-659deae566f6
# ╟─a86b6c8a-9b90-4388-8208-d257ec6951d2
# ╠═d47c1cf4-4844-4c66-8a21-b62ad7740dcf
# ╠═f3f17c65-dab6-46da-83a9-87aad1181524
# ╠═96cf6e3b-1c2b-49f6-979d-0fce28a55a06
# ╠═f7e0f70d-4069-4039-90b6-d3da559bdcb1
# ╠═c7f0582a-8997-4641-b607-0b95a7f73a6c
# ╟─57b35604-727b-4eac-bcbd-ea302e4c79a2
# ╠═3c07c33b-7f91-4e2d-a27a-09a8924f328c
# ╠═b59dd2c0-4c3e-4c3d-a4d4-f8da9dcbe286
# ╠═76e2e659-c708-4a08-8263-def6c3dfa930
# ╟─b62d7318-3a6d-4442-b6cc-3166b87dff8a
# ╠═1bf27aad-b754-4acc-b08d-fe1a86791de8
# ╠═d5beeb75-7ca6-4f28-a63b-6e3070fb4c73
# ╠═5d5f3ca7-260d-4dbb-bde1-59f2bd58c0c9
# ╟─6ac6cce3-09ae-421a-b502-b5c8807bf555
# ╟─cfceae8c-0967-4a47-93f5-2b6e3073624c
# ╠═be090485-7917-464a-8e8c-aec4a2c009d3
# ╠═368f7fc1-0afc-439f-905a-76af04e8ca8d
# ╠═41b0f95e-d806-4167-8199-7392717782f9
# ╠═342cbcff-5a05-4199-b5b5-6cc0fa70c202
# ╠═4249c937-1b33-4aca-8351-d98a19c584fa
# ╠═87df9037-d497-44e8-836d-3c34a58c5d07
# ╟─033473d3-6958-4e7b-8c95-a9e15d882118
# ╟─1ba147d6-2609-4d83-a734-6bf91445a23c
# ╟─7d1f1903-8342-4715-8e88-3e3d063ea5db
# ╟─8790874d-322e-4aea-a60b-4717e653afc9
# ╟─0d122aef-1767-4389-b043-97d07ecef1fe
# ╟─ef346f24-ab46-407b-95a0-2bce10b41ac8
# ╟─c4b39f91-8f0d-4f54-90b4-3c035a5b125f
# ╠═afe1c12f-23d3-47af-bfe8-b59ca7018593
# ╠═666c8d5f-c1c5-4cb4-9364-bf997b18a683
# ╟─b50df30e-f7c3-44f2-b759-3e4ca31003d9
# ╠═3784a5b2-8621-4bac-87a1-4decc25cf4f5
