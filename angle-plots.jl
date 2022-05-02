using GLMakie
import Makie
using LinearAlgebra
using DataFrames
import CSV

df = DataFrame(CSV.File("N:100000,L:10.0,R:106.16233535767984,dc:1.0.csv"))
#df = DataFrame(CSV.File("N:10000,L:10.0,R:49.27619103464619,dc:1.0.csv"))

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
s10 = lines!(ax1,df.links,df.max_compsize_rand)
s11 = lines!(ax1,df.links,df.max_compsize_para)
s12 = lines!(ax1,df.links,df.max_compsize_perp)
s20 = lines!(ax2,df.links,df.global_cluster_coeff_rand)
s21 = lines!(ax2,df.links,df.global_cluster_coeff_para)
s22 = lines!(ax2,df.links,df.global_cluster_coeff_perp)
xlims!.([ax1,ax2], 0, nothing)
ylims!(ax2, 0.0, 0.3)
# ax1.yscale = log10
# ax2.yscale = log10
# ax1.xscale = identity
# ax2.xscale = identity
fig