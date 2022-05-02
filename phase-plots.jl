using GLMakie
import Makie
using LinearAlgebra
using DataFrames
import CSV

df = DataFrame(CSV.File("phase-data.csv"))

# norm so dc_crits = 1.0

Vnormed = @.( df.V / (df.dc_crits^3))
segment_conc_normed = @.(df.N/Vnormed)
segment_length_normed = @.(df.Ls/df.dc_crits)
monomer_conc_normed = @.(segment_conc_normed * segment_length_normed)

fig = Figure()
ax1, s1 = scatter(fig[1,1],segment_length_normed,(segment_conc_normed))
ax2, s2 = scatter(fig[2,1],segment_length_normed,(monomer_conc_normed))
xlims!.([ax1,ax2], 0.0, 10)
ax1.yscale = log10
ax2.yscale = log10
ax1.xscale = identity
ax2.xscale = identity
fig