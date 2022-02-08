module TemplateMatching

using StatsBase
using LinearAlgebra
using OffsetArrays

export crosscorrelate, maxfilter, stack, correlatetemplate, findpeaks, magnitude, findmax_window

include("crosscorrelation.jl")
include("processing.jl")

end