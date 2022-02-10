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

# ╔═╡ 746daaba-22c5-4a47-822a-3a13e687cd37
md"# Data Preparation"

# ╔═╡ cc3e29fc-b063-4421-a2a2-4fe29cb46abe
@bind dirpath TextField(default="../data/2021-01-12_20-25-30")

# ╔═╡ 88625f5f-6e3d-4688-81ae-d6feb3646f8b
samplefreq = 10_000 # in KHz

# ╔═╡ a1e052ff-7660-4cc7-b9f9-9301187a7b3f
readdir(dirpath, join=true)

# ╔═╡ 7485d3e4-5413-4934-85e4-b01ae5596258
listbins(dirpath) = filter(path -> endswith(path, ".bin"), readdir(dirpath, join=true))

# ╔═╡ d07a8c82-7341-4629-b772-b3d8b5324e66
filepaths = listbins(dirpath)

# ╔═╡ 719cdc2e-72db-4790-b606-12c797dbdf76
function convertdata(rawdata)
	data = (convert(Vector{Float32}, ntoh.(reinterpret(Int16, rawdata))))
	data[1:2:end], data[2:2:end]
end

# ╔═╡ 775f53db-5c7c-481a-b1c7-08c6fc49e08d
function readbins(filepaths, nb)
	re = r".+_ch(?P<first_channel_num>[0-9]+)&(?P<second_channel_num>[0-9]+)\.bin"
	data = Dict{Int, Vector{Float32}}()
	for filepath in filepaths
		m = match(re, filepath)
		rawdata = read(filepath, nb)
		first_channel_data, second_channel_data = convertdata(rawdata)
		data[parse(Int, m[:first_channel_num])] = first_channel_data
		data[parse(Int, m[:second_channel_num])] = second_channel_data
	end
	sort(data)
end

# ╔═╡ ba1ac3d7-8376-4955-aeee-1d9daf83964e
md"""

## Bandpass

Enable? $(@bind tofilter CheckBox(default=true))

Lowpass frequency (kHz) $(@bind lopassfreq Slider(10:10:200, default=100, show_value=true))

Highpass frequency (kHz) $(@bind hipassfreq Slider(10:10:800, default=400, show_value=true))

Number of poles in Butterworth filter $(@bind bwpolesnum NumberField(1:10, default=6))
"""

# ╔═╡ 92061d02-57f6-4073-8fea-726a3d82a153
md"""## Resample

Enable? $(@bind toresample CheckBox(default=true))

Decimation factor $(@bind resamplefactor Slider(2:25, default=10, show_value=true))
"""

# ╔═╡ 885e5df3-aa89-4009-88da-23f8e05ca86a
begin
	data = readbins(filepaths, typemax(Int))
	if (tofilter && hipassfreq > lopassfreq)
		responsetype = Bandpass(lopassfreq, hipassfreq, fs=samplefreq)
		designmethod = Butterworth(bwpolesnum)
		Threads.@threads for n in collect(keys(data))
			data[n] = filtfilt(digitalfilter(responsetype, designmethod), data[n])
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
Plot scale $(@bind plotnpts Slider(100:100:5000, default=1000, show_value=true))

Plot window position $(@bind plotpos Slider(0:0.01:100, default=50, show_value=true))
"""

# ╔═╡ c717e19e-8007-4f45-ada5-e2879b238e0d
function plotstream(data, npts, center)
	plots = []
	numchannels = length(data)
	for (nch, series) in data
		a = max(first(axes(series, 1)), center - div(npts, 2))
		b = min(last(axes(series, 1)), center + div(npts, 2))
		ylim = extrema(series)
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

# ╔═╡ dd071ca8-ff4e-4af2-9d44-48e7947c3431
let center = round(Int, (plotpos / 100) * minimum(size(series, 1) for (_, series) in data)) 
	plotstream(data, plotnpts, center)
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
# ╠═88625f5f-6e3d-4688-81ae-d6feb3646f8b
# ╠═a1e052ff-7660-4cc7-b9f9-9301187a7b3f
# ╠═7485d3e4-5413-4934-85e4-b01ae5596258
# ╠═d07a8c82-7341-4629-b772-b3d8b5324e66
# ╠═719cdc2e-72db-4790-b606-12c797dbdf76
# ╠═775f53db-5c7c-481a-b1c7-08c6fc49e08d
# ╟─ba1ac3d7-8376-4955-aeee-1d9daf83964e
# ╟─92061d02-57f6-4073-8fea-726a3d82a153
# ╟─885e5df3-aa89-4009-88da-23f8e05ca86a
# ╟─67b94e20-6332-4993-bdef-8cc29e51b1ea
# ╟─c717e19e-8007-4f45-ada5-e2879b238e0d
# ╟─dd071ca8-ff4e-4af2-9d44-48e7947c3431
# ╟─60930b7e-c200-4219-929b-58e33cd8b939
# ╟─a4728746-a8a0-4bbb-9e9e-228eee132fdf
# ╟─246153f0-c3ed-4cec-b88c-3d8bb5790c17
# ╟─db92c8ac-9fce-4e36-bb2f-de917654af00
# ╠═01c58f23-0da5-4d5f-9530-c96d43c655a0
# ╟─a86b6c8a-9b90-4388-8208-d257ec6951d2
# ╠═908feaef-7304-48b4-b8a5-bea3ab7d8e43
