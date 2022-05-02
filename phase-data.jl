using LinearAlgebra
using DataFrames
import CSV
using ProgressMeter

include("util.jl")

trials = 1000

N = 1_000_000

V = 4.0*N

R = cbrt(3//4*V/pi)

Ls    = LinRange(0.1,4.0,trials)
dcmax = 1.5

dc_crits = fill(NaN, trials)

p = Progress(trials)

Threads.@threads for i in 1:trials
    L = Ls[i]
    dc_crits[i] = critical_pt(N,L,R,dcmax)
    next!(p)
end

# save data
df = DataFrame(;
    N,
    V,
    Ls,
    dcmax,
    dc_crits,
)

CSV.write("phase-data.csv",df)


