module TemplateMatching

using StatsBase
using LinearAlgebra
using OffsetArrays

export crosscorrelate, maxfilter, stack, correlatetemplate, findpeaks, magnitude

include("crosscorrelation.jl")
include("processing.jl")

end