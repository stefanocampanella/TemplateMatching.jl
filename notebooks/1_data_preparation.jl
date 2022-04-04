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

# ╔═╡ e81954be-8854-46eb-9f26-0e2d5c8896e5
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()
end

# ╔═╡ 8b4426b8-588b-49ad-83c5-a79ed698d704
using PlutoUI

# ╔═╡ 18fbfdc6-9f0b-451f-a6df-fb602725a90c
using DSP

# ╔═╡ 8382586d-170f-45fd-82eb-e1b526b2cf1c
using JLD2

# ╔═╡ f07ab63f-dfb4-493c-a40c-6b3363a858c3
using Plots

# ╔═╡ 58ba75d2-ab31-43b5-a266-df944fd79530
using StatsBase

# ╔═╡ 746daaba-22c5-4a47-822a-3a13e687cd37
md"# Data Preparation"

# ╔═╡ cc3e29fc-b063-4421-a2a2-4fe29cb46abe
@bind dirpath TextField(default="../data/2021-01-12_20-25-30")

# ╔═╡ d07a8c82-7341-4629-b772-b3d8b5324e66
filepaths = filter(path -> endswith(path, ".bin"), readdir(dirpath, join=true))

# ╔═╡ 719cdc2e-72db-4790-b606-12c797dbdf76
function convertdata(rawdata, eltype::Type{T}) where T <: AbstractFloat
	data = (convert(Vector{eltype}, ntoh.(reinterpret(Int16, rawdata))))
	data[1:2:end], data[2:2:end]
end

# ╔═╡ 775f53db-5c7c-481a-b1c7-08c6fc49e08d
function readbins(filepaths, nb, eltype::Type{T}) where T <: AbstractFloat
	re = r".+_ch(?P<first_channel_num>[0-9]+)&(?P<second_channel_num>[0-9]+)\.bin"
	data = Dict{Int, Vector{eltype}}()
	for filepath in filepaths
		m = match(re, filepath)
		rawdata = read(filepath, nb)
		first_channel_data, second_channel_data = convertdata(rawdata, eltype)
		data[parse(Int, m[:first_channel_num])] = first_channel_data
		data[parse(Int, m[:second_channel_num])] = second_channel_data
	end
	sort(data)
end

# ╔═╡ ba1ac3d7-8376-4955-aeee-1d9daf83964e
md"""

## Bandpass

Enable? $(@bind tofilter CheckBox(default=true))

Lowpass frequency (kHz) $(@bind lopassfreq Slider(1:2000, default=50, show_value=true))

Highpass frequency (kHz) $(@bind hipassfreq Slider(1:2000, default=600, show_value=true))

Number of poles in Butterworth filter $(@bind bwpolesnum NumberField(1:10, default=4))
"""

# ╔═╡ 92061d02-57f6-4073-8fea-726a3d82a153
md"""## Resample

Enable? $(@bind toresample CheckBox(default=true))

Decimation factor $(@bind resamplefactor Slider(2:25, default=5, show_value=true))
"""

# ╔═╡ 88625f5f-6e3d-4688-81ae-d6feb3646f8b
samplefreq = 10_000 # in KHz

# ╔═╡ 885e5df3-aa89-4009-88da-23f8e05ca86a
begin
	data = readbins(filepaths, typemax(Int), Float64)
	if (tofilter && hipassfreq > lopassfreq)
		responsetype = Bandpass(lopassfreq, hipassfreq, fs=samplefreq)
		designmethod = Butterworth(bwpolesnum)
		Threads.@threads for n in collect(keys(data))
			ys = data[n]
			xs = axes(data[n], 1)
			β = (mean(xs .* ys) - mean(xs) * mean(ys)) / std(xs, corrected=false)
			α = mean(ys) - β * mean(xs)
			δs = @. ys - α - β * xs
			data[n] = filtfilt(digitalfilter(responsetype, designmethod), δs)
		end
	end
	if toresample
		Threads.@threads for n in collect(keys(data))
			data[n] = resample(data[n], 1 //resamplefactor)
		end
	end
	data
end

# ╔═╡ 67b94e20-6332-4993-bdef-8cc29e51b1ea
md"""
Plot x-range $(@bind npts Slider(100:100:5000, default=1000, show_value=true))

Plot window position $(@bind center_percent Slider(0:0.01:100, default=50, show_value=true))
"""

# ╔═╡ dd071ca8-ff4e-4af2-9d44-48e7947c3431
let 
	center = round(Int, (center_percent / 100) * minimum(size(series, 1) for (_, series) in data)) 
	plots = []
	numchannels = length(data)
	for (nch, series) in data
		a = max(firstindex(series, 1), center - div(npts, 2))
		b = min(lastindex(series, 1), center + div(npts, 2))
		ylim = extrema(view(series, a:b))
		push!(plots, plot(a:b, view(series, a:b),
			  ylim=ylim,
			  title="Channel $nch",
			  titlealign=:left,
			  showaxis=false,
			  ticks=nothing,
			  label=nothing))
	end
	plot(plots..., layout=(numchannels, 1), size=(800, 100numchannels))
end

# ╔═╡ 60930b7e-c200-4219-929b-58e33cd8b939
new_samplefreq = samplefreq // resamplefactor

# ╔═╡ a4728746-a8a0-4bbb-9e9e-228eee132fdf
pgrams = Dict(channel => periodogram(series, fs=new_samplefreq) for (channel, series) in data)

# ╔═╡ 246153f0-c3ed-4cec-b88c-3d8bb5790c17
md"""
Which channel? $(@bind pgramplotch Select(sort(collect(keys(pgrams)))))

Downsampling $(@bind pgram_downsample Slider(10:10:100, default=50, show_value=true))
"""

# ╔═╡ db92c8ac-9fce-4e36-bb2f-de917654af00
plot(resample(pgrams[pgramplotch].freq, 1//pgram_downsample),
	 resample(pgrams[pgramplotch].power, 1//pgram_downsample), 
     title="Channel $pgramplotch", label=nothing, dpi=300)

# ╔═╡ 01c58f23-0da5-4d5f-9530-c96d43c655a0
bad_channels = [1, 2, 3]

# ╔═╡ a86b6c8a-9b90-4388-8208-d257ec6951d2
@bind outputpath TextField(default="../data/2021-01-12_20-25-30/data.jld2")

# ╔═╡ 908feaef-7304-48b4-b8a5-bea3ab7d8e43
save(outputpath, "data", filter(!(x -> x.first in bad_channels), data))

# ╔═╡ Cell order:
# ╟─746daaba-22c5-4a47-822a-3a13e687cd37
# ╠═e81954be-8854-46eb-9f26-0e2d5c8896e5
# ╠═8b4426b8-588b-49ad-83c5-a79ed698d704
# ╠═18fbfdc6-9f0b-451f-a6df-fb602725a90c
# ╠═8382586d-170f-45fd-82eb-e1b526b2cf1c
# ╠═f07ab63f-dfb4-493c-a40c-6b3363a858c3
# ╠═cc3e29fc-b063-4421-a2a2-4fe29cb46abe
# ╠═d07a8c82-7341-4629-b772-b3d8b5324e66
# ╠═719cdc2e-72db-4790-b606-12c797dbdf76
# ╠═775f53db-5c7c-481a-b1c7-08c6fc49e08d
# ╟─ba1ac3d7-8376-4955-aeee-1d9daf83964e
# ╠═92061d02-57f6-4073-8fea-726a3d82a153
# ╠═88625f5f-6e3d-4688-81ae-d6feb3646f8b
# ╠═58ba75d2-ab31-43b5-a266-df944fd79530
# ╠═885e5df3-aa89-4009-88da-23f8e05ca86a
# ╟─67b94e20-6332-4993-bdef-8cc29e51b1ea
# ╟─dd071ca8-ff4e-4af2-9d44-48e7947c3431
# ╟─60930b7e-c200-4219-929b-58e33cd8b939
# ╟─a4728746-a8a0-4bbb-9e9e-228eee132fdf
# ╟─246153f0-c3ed-4cec-b88c-3d8bb5790c17
# ╟─db92c8ac-9fce-4e36-bb2f-de917654af00
# ╠═01c58f23-0da5-4d5f-9530-c96d43c655a0
# ╟─a86b6c8a-9b90-4388-8208-d257ec6951d2
# ╠═908feaef-7304-48b4-b8a5-bea3ab7d8e43
