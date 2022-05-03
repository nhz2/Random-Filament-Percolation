using CairoMakie
import Makie
using LinearAlgebra
using DataFrames
import CSV

df = DataFrame(CSV.File("N100000 L10 R106 dc1.csv"))
#df = DataFrame(CSV.File("N:10000,L:10.0,R:49.27619103464619,dc:1.0.csv"))
N = 100000

begin
fig = Figure( fontsize = 20)
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
s10 = lines!(ax1,df.links,df.max_compsize_rand./N; label="Random")
s11 = lines!(ax1,df.links,df.max_compsize_para./N; label="∥ first")
s12 = lines!(ax1,df.links,df.max_compsize_perp./N; label="⟂ first")
s20 = lines!(ax2,df.links,df.global_cluster_coeff_rand; label="Random")
s21 = lines!(ax2,df.links,df.global_cluster_coeff_para; label="∥ first")
s22 = lines!(ax2,df.links,df.global_cluster_coeff_perp; label="⟂ first")
xlims!.([ax1,ax2], 0, nothing)
ax2.xlabel = "Links added"
ylims!(ax1, 0.0, 1)
ax1.ylabel = "Max component size / N"
ylims!(ax2, 0.0, 0.3)
ax2.ylabel = "Global clustering coeff"
ax1.yticks = LinearTicks(4)
ax2.yticks = LinearTicks(3)
axislegend(ax1; position=:rb)
#Legend(fig[1:2, 2],[s11,s10,s12],["∥ first","Random","⟂ first"])
fig
end

save("angle-plots.svg",fig)