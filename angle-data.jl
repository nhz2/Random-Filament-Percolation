using LinearAlgebra
using DataFrames
import CSV
using ProgressMeter

include("util.jl")

N = 100_000
V = 10^1.7*N
R = cbrt(3//4*V/pi)
L = 10.0
distancecutoff = 1.0
@time df=angle_based_bonding(N,L,R,distancecutoff,10)

CSV.write("N:$N,L:$L,R:$R,dc:$distancecutoff.csv",df)