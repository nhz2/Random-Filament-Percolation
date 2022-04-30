# Utility functions
# Nathan Zimmerberg 30 APR 2022

using StaticArrays
using LinearAlgebra
import CellListMap
using Random

"""
Return a 2 by N array of 3D points, each column represents a line segment.
All line segments have length L, and are randomly placed in fully in a sphere of radius R.
"""
function generate_segments(N, R, L)::Array{SVector{3, Float64}, 2}
    a = Array{SVector{3, Float64}, 2}(undef,2,N)
    L < 2R || error("line segments cannot fit in the sphere")
    i = 0
    Rsqr = R^2
    halfL = (1//2)*L
    while i < N
        #try and add a line segment
        center = 2R*(rand(SVector{3, Float64}) .- 0.5)
        dir = normalize(randn(SVector{3, Float64}))
        point1 = center + halfL*dir
        point2 = center - halfL*dir
        if (point1 ⋅ point1) < Rsqr && (point2 ⋅ point2) < Rsqr
            i += 1
            a[1,i] = point1
            a[2,i] = point2
        end
    end
    return a
end


"""
    linesegment_linesegment_dist2(P0,P1,Q0,Q1)

Return the squared distance between two line segments.
Using the simple algorithm and some comments from
https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2022
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 6.0.2022.01.06
Ignores the case of degenerate line segments.
"""
function linesegment_linesegment_dist2(P0,P1,Q0,Q1)
    P = P1-P0
    Q = Q1-Q0
    P0mQ0 = P0 - Q0
    a = P ⋅ P
    b = P ⋅ Q
    c = Q ⋅ Q
    d = P ⋅ P0mQ0
    e = Q ⋅ P0mQ0
    f = P0mQ0 ⋅ P0mQ0
    Δ = a*c - b^2
    #assuming both segments are not zero length
    #critical points
    s0 = clamp(-d/a,0,1)
    s1 = clamp((b-d)/a,0,1)
    t0 = clamp(e/c,0,1)
    t1 = clamp((b+e)/c,0,1)
    sbar = clamp((b*e - c*d)/(Δ+eps(one(Δ))),0,1)
    tbar = clamp((a*e - b*d)/(Δ+eps(one(Δ))),0,1)
    r(s,t) = a*s^2 - 2b*s*t + c*t^2 + 2d*s - 2e*t
    return max(0,f+min(
        r(s0,0),
        r(s1,1),
        r(0,t0),
        r(1,t1),
        r(sbar,tbar),
    ))
end


struct LinkableSegments
    id1::Int32
    id2::Int32
    d2::Float32 #squared distance
    cosθ::Float32 #cos angle 
end

"""
Return a vector of LinkableSegments that are within distancecutoff
"""
function generate_neighborlist(linesegments::Array{SVector{3, Float64}, 2}, distancecutoff)::Vector{LinkableSegments}
    centers = (linesegments[1,:] + linesegments[2,:])/2
    maxlength = maximum(norm,(linesegments[1,:] - linesegments[2,:]))
    centercutoff = distancecutoff + maxlength
    box = CellListMap.Box(CellListMap.limits(centers),centercutoff)
    cl = CellListMap.CellList(centers,box)
    nl = Vector{LinkableSegments}()
    CellListMap.map_pairwise_serial!(nl,box,cl) do x,y,i,j,d2,output
        i0 = linesegments[1,i]
        i1 = linesegments[2,i]
        j0 = linesegments[1,j]
        j1 = linesegments[2,j]
        d2 = linesegment_linesegment_dist2(i0,i1,j0,j1)
        if d2 ≤ distancecutoff^2
            i_dir = normalize(i1-i0)
            j_dir = normalize(j1-j0)
            cosθ = i_dir ⋅ j_dir
            push!(output,LinkableSegments(i,j,d2,cosθ))
        end
        output
    end
    nl
end


"""
Error if a neighbor list is incorrect.
"""
function check_neighborlist(nl::Vector{LinkableSegments}, linesegments::Array{SVector{3, Float64}, 2}, distancecutoff)
    neighbor_set = Set(nl)
    N = size(linesegments,2)
    for i in 1:N-1
        for j in i+1:N
            i0 = linesegments[1,i]
            i1 = linesegments[2,i]
            j0 = linesegments[1,j]
            j1 = linesegments[2,j]
            d2 = linesegment_linesegment_dist2(i0,i1,j0,j1)
            if d2 < distancecutoff^2
                i_dir = normalize(i1-i0)
                j_dir = normalize(j1-j0)
                cosθ = i_dir ⋅ j_dir
                i2j = LinkableSegments(i,j,d2,cosθ)
                j2i = LinkableSegments(j,i,d2,cosθ)
                (i2j in neighbor_set) || (j2i in neighbor_set) || error("$i, $j, not in neighbor list")
            end
        end
    end
end