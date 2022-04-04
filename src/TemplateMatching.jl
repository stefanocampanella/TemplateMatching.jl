module TemplateMatching

using StatsBase
using LinearAlgebra
using OffsetArrays
using Interpolations
using Optim

export crosscorrelate, maxfilter, stack, correlatetemplate, findpeaks, magnitude, findmax_window, estimatetoa

include("crosscorrelation.jl")
include("processing.jl")

end