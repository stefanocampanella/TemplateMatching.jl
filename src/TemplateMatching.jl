module TemplateMatching

using StatsBase
using Primes
using FFTW
using LinearAlgebra
using OffsetArrays
using Interpolations
using Optim

export crosscorrelate, maxfilter, stack, correlatetemplate, findpeaks, 
       magnitude, findmax_window, subsampleshift, estimatetoa, residue_rms, 
       locate

include("crosscorrelation.jl")
include("processing.jl")
include("utils.jl")

end