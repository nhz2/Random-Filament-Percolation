# Utility functions
# Nathan Zimmerberg 30 APR 2022

using StaticArrays
using LinearAlgebra
import CellListMap
using Random
using Graphs
using DataFrames

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

"""
Percolation Algorithm from Lecture041922
"""
struct QuickUnion
    parents::Vector{Int32}
    sizes::Vector{Int32}
    pathbuffer::Vector{Int32}
end

function QuickUnion(n::Integer)
    pathbuffer = Vector{Int32}()
    sizehint!(pathbuffer,n)
    QuickUnion(zeros(n),ones(n),pathbuffer)
end

function findroot(qu::QuickUnion,i::Int32)::Int32
    parents = qu.parents
    pathbuffer = qu.pathbuffer
    parenti = parents[i]
    if !iszero(parenti)
        while true
            nextparent = parents[parenti]
            if !iszero(nextparent)
                push!(pathbuffer, i)
                i = parenti
                parenti = nextparent
            else
                # parenti is the root
                i = parenti
                break
            end
        end
    end
    root = i
    # i is now root
    # go through and set all parents to root
    for pathi in pathbuffer
        parents[pathi] = root
    end
    empty!(pathbuffer)
    return root
end

"""
Return the size of the merged clusters, or 0 if no clusters were merged
"""
function connect!(qu::QuickUnion,a::Int32,b::Int32)::Int32
    roota = findroot(qu,a)
    rootb = findroot(qu,b)
    if roota == rootb
        return 0
    end
    #add pointer from smaller cluster to larger
    parents = qu.parents
    sizes = qu.sizes
    sizea = sizes[roota]
    sizeb = sizes[rootb]
    if sizea < sizeb
        parents[roota] = rootb
        sizes[roota] = 0
        sizes[rootb] = sizea + sizeb
    else
        parents[rootb] = roota
        sizes[rootb] = 0
        sizes[roota] = sizea + sizeb
    end
    return sizea + sizeb
end



"""
Return the critical distancecutoff 
where 10% of line segments are in the largest connected component, 
or nothing if the largest connected component is never large enough
"""
function critical_pt(N, L, R, maxdistancecutoff)::Float64
    segs = generate_segments(N, R, L)
    nl = generate_neighborlist(segs, maxdistancecutoff)
    #@show length(nl)
    if length(nl)+1 < 0.1*N
        return NaN
    end
    # sort nl from closest to farthest distance
    sort!(nl; by=(x->x.d2))
    max_cluster_size = 1
    # create quick union struct
    qu = QuickUnion(N)
    for neighbor::LinkableSegments in nl
        newsize = connect!(qu,neighbor.id1,neighbor.id2)
        max_cluster_size = max(max_cluster_size,newsize)
        if max_cluster_size ≥ 0.1*N
            return √(neighbor.d2)
        end
    end
    #@show max_cluster_size
    return NaN
end



"""
Return a DataFrame with columns of links, compsize_para, compsize_perp, global_cluster_coeff_para, global_cluster_coeff_perp

First create random line segments and a neighbor list
Add links using two strategies.
    1. para: add pairs of parallel or anti parallel line segments first.
    2. perp: add pairs of perpendicular line segments first.
The angle is determined by the absolute value of the dot product of the two line segment direction vectors
"""
function angle_based_bonding(N, L, R, distancecutoff)
    segs = generate_segments(N, R, L)
    nl = generate_neighborlist(segs, maxdistancecutoff)
    error("todo")
end