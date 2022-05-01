import GeometryBasics
using GLMakie
import Makie
using LinearAlgebra

include("util.jl")

"""
Return a GeometryBasics Mesh of a capsule shape, centered at the origin, 
length l, radius r, aligned in the z direction.
Using examples from Makie.jl
https://makie.juliaplots.org/stable/examples/plotting_functions/mesh/
"""
function capsule_mesh(l,r,n)
    φs = [LinRange(pi, pi/2, n); LinRange(pi/2, 0, n)]
    zoffsets = [fill(-l/2, n); fill(+l/2, n)]
    θs = LinRange(0, 2pi, 4*n)
    
    points = Vector{Point3f}()
    normals = Vector{Point3f}()
    for (φ, zoffset) in zip(φs,zoffsets)
        for θ in θs
            x = sin(φ)*cos(θ)
            y = sin(φ)*sin(θ)
            z = cos(φ)
            push!(normals,Vec3f(x,y,z))
            push!(points,Point3f(r*x,r*y,r*z + zoffset))
        end
    end

    # The coordinates form a matrix, so to connect neighboring vertices with a face
    # we can just use the faces of a rectangle with the same dimension as the matrix:
    faces = GeometryBasics.decompose(GeometryBasics.TriangleFace{GeometryBasics.GLIndex}, GeometryBasics.Tesselation(GeometryBasics.Rect(0, 0, 1, 1), (4n,2n)))
    # @show faces
    # Normals of a centered sphere are easy, they're just the vertices normalized.
    GeometryBasics.Mesh(GeometryBasics.meta(points; normals), faces)
end
N = 10000
gb_mesh = capsule_mesh(10.0,0.5,10)

linesegments = generate_segments(N,100.0,10.0)

centers = (linesegments[1,:] + linesegments[2,:])/2

zaxis = normalize.(linesegments[1,:] - linesegments[2,:])

colors = fill(:blue,N)

scene = Scene(backgroundcolor=:white)
cam = cam3d!(scene)

# w, h = size(scene)
# nearplane = 0.1f0
# farplane = 100f0
# aspect = Float32(w / h)
# cam.projection[] = Makie.perspectiveprojection(
#     45f0, 
#     aspect, 
#     nearplane, 
#     farplane
# )


meshscatter!(scene, centers;
    marker = gb_mesh,
    markersize = 1.0,
    rotation = zaxis,
    color = colors, 
    # diffuse = Vec3f(0.3), 
    # specular = Vec3f(0.2),
    # shininess = 32.0,
    ssao = true,
)

center!(scene)

scene
# f, ax, pl = mesh(gb_mesh;
#     color = :blue, 
#     # diffuse = Vec3f(0.1), 
#     # specular = Vec3f(0.1),
#     # shininess = 32.0,
#     ssao = true,
#     show_axis = false,
# )
# wireframe!(ax, gb_mesh, color=(:black, 0.2), linewidth=2, transparency=true)
# f