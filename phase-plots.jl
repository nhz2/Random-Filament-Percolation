using CairoMakie
import Makie
using LinearAlgebra
using DataFrames
import CSV

df = DataFrame(CSV.File("phase-data.csv"))

fig = Figure( fontsize = 25)
ax1, s1 = scatter(fig[1,1],df.Ls,df.dc_crits;
    markersize=3,
)
xlims!(ax1, 0, 4)
ax1.yscale = identity
ax1.xscale = identity
ax1.xlabel = "L"
ax1.ylabel = "Critical cutoff distance"
ax1.title = "Critical cutoff distance vs L at N = 10⁶ and V = 4N"
fig
save("Critical-points-raw.svg",fig)

# norm so dc_crits = 1.0

Vnormed = @.( df.V / (df.dc_crits^3))
segment_conc_normed = @.(df.N/Vnormed)
segment_length_normed = @.(df.Ls/df.dc_crits)
monomer_conc_normed = @.(segment_conc_normed * segment_length_normed)

begin
fig = Figure( fontsize = 16)
ax1, s1 = scatter(fig[1,1],segment_length_normed,segment_conc_normed;
    markersize=3,
)
ax2, s2 = scatter(fig[2,1],segment_length_normed,monomer_conc_normed;
    markersize=3,
)
xlims!.([ax1,ax2], 0.0, 10)
ax2.xlabel = "Aspect ratio: L/dc"
ax1.yscale = log10
ax2.yscale = log10
ax1.xscale = identity
ax2.xscale = identity
ax1.ylabel = "Segment conc: N / (V/dc³)"
ax2.ylabel = "Monomer conc: N*(L/dc) / (V/dc³)"
ax1.xticks = LinearTicks(10)
ax2.xticks = LinearTicks(10)
fig
end

save("Critical-points-normed.svg",fig)