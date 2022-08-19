module TemplateMatching

using StatsBase
using Primes
using FFTW
using CUDA
using LinearAlgebra
using OffsetArrays
using Interpolations

export crosscorrelate, maxfilter, stack, correlatetemplate, findpeaks, 
       relative_magnitude, findmax_window, subsampleshift, estimatetoa

include("crosscorrelation.jl")
include("utils.jl")

end