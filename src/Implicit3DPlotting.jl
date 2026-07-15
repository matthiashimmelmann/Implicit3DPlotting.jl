module Implicit3DPlotting

#=
    This is a package for plotting implicitly defined space curves and surfaces. It plots the objects in a single window
    `GLMakie.Figure` and `WGLMakie.Scene` are used for this. If multiple simultaneous windows are necessary, we recommend
    either using the web-based openGL backend in Jupyter notebooks. If that does not work, you can ask for this feature; 
    it did not seem necessary up until now.

    author:     Matthias Himmelmann
    email:      matthias.himmelmann (at) outlook.de
    website:    matthiashimmelmann.github.io/
=#

import GLMakie: Figure, wireframe!, linesegments!, mesh!, Point3f, DirectionalLight, LScene, RGBf, Axis3, hidedecorations!, hidespines!
import GLMakie.GeometryBasics: Mesh, Vec3f, decompose, TriangleFace, Point, normal_mesh
import Meshing: MarchingCubes, MarchingTetrahedra, isosurface
import LinearAlgebra: pinv, norm, nullspace, cross, dot
import HomotopyContinuation: evaluate, differentiate, @var, Expression
import IterTools: product

export plot_implicit_surface,
       plot_implicit_surface!,
       plot_implicit_curve,
       plot_implicit_curve!
#INFO: The following packages are not maintained by me. Find them here: https://github.com/MakieOrg/Makie.jl


function plot_implicit_surface!(
    ax::Axis3,
    f::Function;
    xlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (-3,3),
    ylims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (-3,3),
    zlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (-3,3),
    color::Union{Symbol,Tuple{Symbol,Float64}} = :steelblue,
    transparency::Bool = false,
    samples::Tuple{Int,Int,Int}=(35,35,35),
    wireframe::Bool=false,
    MarchingModeIsCubes::Bool=true,
    color_gradient=:turbo,
    color_mapping=nothing,
    cutoffmap=nothing,
    diffuse = Vec3f(0.8),
    specular = Vec3f(1.1),
    shininess = 60f0,
    backlight = 5f0,
    kwargs...
)

    try
        f([1,1,1])
    catch
        try
            f(1,1,1)
            f = x->f(x[1],x[2],x[3])
        catch
            throw("The specified function `f` does not take 3 arguments.")
        end
    end
    algo = MarchingModeIsCubes ? MarchingCubes(; iso=0) : MarchingTetrahedra(; iso=0)
    ξx=xlims[1]:(xlims[2]-xlims[1])/samples[1]:xlims[2]
    ξy=ylims[1]:(ylims[2]-ylims[1])/samples[2]:ylims[2]
    ξz=zlims[1]:(zlims[2]-zlims[1])/samples[3]:zlims[2]
    f_values = [f([x,y,z]) for x in ξx, y in ξy, z in ξz]
    vertices, facets = isosurface(f_values, algo, ξx, ξy, ξz)
    implicit_mesh = Mesh(Point3f.(vertices), TriangleFace.(facets))

    if wireframe
        wireframe!(ax, implicit_mesh; 
            color=color, 
            transparency=transparency, 
            kwargs...
        )
    else
        if isnothing(color_mapping)
            mesh!(ax, normal_mesh(implicit_mesh); 
                color = color, 
                diffuse = diffuse,
                specular = specular,
                shininess = shininess,
                backlight = backlight,
                transparency = transparency,
                kwargs...
            )
        else
            colors = [color_mapping(v) for v in vertices]
            colors = (colors .- minimum(colors)) ./ (maximum(colors)-minimum(colors))
            mesh!(ax, normal_mesh(implicit_mesh);
                color = colors, 
                diffuse = diffuse,
                specular = specular,
                shininess = shininess,
                backlight = backlight,
                transparency=transparency,
                colormap=color_gradient, 
                kwargs...
            )
        end        
    end
    return nothing
end

#=
Plots an implicitly defined surface.
=#
function plot_implicit_surface(
    f::Function;
    x_min = -3.0,
    x_max = 3.0,
    y_min = x_min,
    y_max = x_max,
    z_min = x_min,
    z_max = x_max,
    xlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (x_min, x_max),
    ylims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (y_min, y_max),
    zlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (z_min, z_max),
    show_axis = false,
    resolution=(800,800),
    aspect=(1.,1.,1.),
    samples::Tuple{Int,Int,Int}=(35,35,35),
    transparency=true,
    fontsize=17,
    diffuse = Vec3f(1.1),
    specular = Vec3f(0.8),
    shininess = 60f0,
    backlight = 4f0,
    kwargs...
)
    global fig = Figure(size=resolution)
    ax = Axis3(fig[1,1], aspect=aspect)
    if !show_axis
        hidedecorations!(ax)
        hidespines!(ax)
    end
    plot_implicit_surface!(ax, f; 
        xlims=xlims, 
        ylims=ylims, 
        zlims=zlims, 
        transparency=transparency, 
        diffuse=diffuse,
        samples=samples,
        specular=specular, 
        shininess=shininess, 
        backlight=backlight, 
        kwargs...
    )
    return(fig)
end

"""
Adds a space curve, implicitly defined by two equations, to a scene.
"""
function plot_implicit_curve!(
    ax,
    f::Function,
    g::Function;
    color::Union{Symbol,Tuple{Symbol,Float64}} = :steelblue,
    linewidth::Union{Float64,Int}=2.25,
    MarchingModeIsCubes::Bool = true,
    samples::Tuple{Int,Int,Int}=(10,10,10),
    step_size=0.025,
    transparency::Bool = true,
    cutoffmap=nothing,
    miniter=5,
    xlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (-2.,2),
    ylims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (-2.,2),
    zlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (-2.,2),
    kwargs...
)
    try
        f([1,1,1])
    catch
        try
            f(1,1,1)
            f = x->f(x[1],x[2],x[3])
        catch
            throw("The specified function `f` does not take 3 arguments.")
        end
    end
     try
        g([1,1,1])
    catch
        try
            g(1,1,1)
            g = x->g(x[1],x[2],x[3])
        catch
            throw("The specified function `g` does not take 3 arguments.")
        end
    end

    lines=[]

    try
        @var x[1:3]
        poly_f = Expression(f(x))
        poly_g = Expression(g(x))
        point_samples, lines = [], []
        poly_sys = [poly_f, poly_g]
        jac = differentiate(poly_sys, x)
        for i in 1:Int(round(sum(samples)/3))
            lin_space = (rand(Float64,3) .- 0.5)'*x
            lin_jac = differentiate(vcat(poly_sys, lin_space), x)
            for _ in 1:4
                start_point = (i<=2 ? 0.025 : 0.5)*[xlims[2]-xlims[1], ylims[2]-ylims[1], zlims[2]-zlims[1]] .* (rand(Float64,3) .- 0.5)
                try
                    q = newtoncorrect(vcat(poly_sys, lin_space), x, lin_jac, start_point)
                    if !any(t->isapprox(norm(t-q),0; atol=1e-10), point_samples)
                        push!(point_samples, q)
                    end
                catch
                    continue
                end
            end
        end

        for point_sample in point_samples
            if any(tl->isapprox(tl[1], point_sample, atol=step_size) || isapprox(tl[2], point_sample, atol=step_size), lines)
                continue
            end
            one_line = [Point3f(point_sample)]
            q = Base.copy(point_sample)
            prev_flex = nullspace(evaluate(jac, x=>point_sample))[:,1]
            i = 1
            while true
                q, prev_flex = euler_step(jac, x, step_size, prev_flex, q)
                q = newtoncorrect(poly_sys, x, jac, q)
                push!(one_line, Point3f(q))
                i=i+1
                if i<miniter
                    continue
                end
                if (cutoffmap!=nothing && !cutoffmap(q)) || (q[1]<xlims[1] || q[1]>xlims[2] || q[2]<ylims[1] || q[2]>ylims[2] || q[3]<zlims[1] || q[3]>zlims[2] ) || isapprox(norm(q-point_sample),0; atol=step_size)
                    break
                end
            end

            for i in 1:length(one_line)-1
                if !any(tl->isapprox(tl[1], one_line[i], atol=step_size/2) && isapprox(tl[2], one_line[i+1], atol=step_size/2), lines) && !any(tl->isapprox(tl[2], one_line[i], atol=step_size/2) && isapprox(tl[1], one_line[i+1], atol=step_size/2), lines)
                    push!(lines, [one_line[i], one_line[i+1]])
                end
            end

            if isapprox(point_sample, one_line[end], atol=step_size)
                continue
            end

            # Also run the algorithm in the opposite direction
            one_line = [Point3f(point_sample)]
            q = Base.copy(point_sample)
            prev_flex = -nullspace(evaluate(jac, x=>point_sample))[:,1]
            i = 1
            while true
                q, prev_flex = euler_step(jac, x, step_size, prev_flex, q)
                q = newtoncorrect(poly_sys, x, jac, q)
                push!(one_line, Point3f(q))
                i=i+1
                if i<miniter
                    continue
                end
                if (cutoffmap!=nothing && !cutoffmap(q)) || (q[1]<xlims[1] || q[1]>xlims[2] || q[2]<ylims[1] || q[2]>ylims[2] || q[3]<zlims[1] || q[3]>zlims[2] ) || isapprox(norm(q-point_sample),0; atol=step_size)
                    break
                end
            end

            for i in 1:length(one_line)-1
                if !any(tl->isapprox(tl[1], one_line[i], atol=step_size/2) && isapprox(tl[2], one_line[i+1], atol=step_size/2), lines) && !any(tl->isapprox(tl[2], one_line[i], atol=step_size/2) && isapprox(tl[1], one_line[i+1], atol=step_size/2), lines)
                    push!(lines, [one_line[i], one_line[i+1]])
                end
            end

        end

    catch e
        @warn "The given functions `f` and `g` could not be cast to polynomials. Error code: $(e)."        

        algo = MarchingModeIsCubes ? MarchingCubes(; iso=0) : MarchingTetrahedra(; iso=0)
        ξx=xlims[1]:(xlims[2]-xlims[1])/samples[1]:xlims[2]
        ξy=ylims[1]:(ylims[2]-ylims[1])/samples[2]:ylims[2]
        ξz=zlims[1]:(zlims[2]-zlims[1])/samples[3]:zlims[2]
        f_values = [f([x,y,z]) for x in ξx, y in ξy, z in ξz]
        f_vertices, f_facets = isosurface(f_values, algo, ξx, ξy, ξz)
        g_values = [g([x,y,z]) for x in ξx, y in ξy, z in ξz]
        g_vertices, g_facets = isosurface(g_values, algo, ξx, ξy, ξz)

        for triangle_f in f_facets, triangle_g in g_facets
            vertices_for_f = [collect(f_vertices[tri]) for tri in triangle_f]
            vertices_for_g = [collect(g_vertices[tri]) for tri in triangle_g]
            if (min(sum((vertices_for_f[1]-vertices_for_g[1]).^2), sum((vertices_for_f[2]-vertices_for_g[2]).^2), sum((vertices_for_f[3]-vertices_for_g[3]).^2))
                    < max(sum((vertices_for_f[1]-vertices_for_f[2]).^2), sum((vertices_for_f[2]-vertices_for_f[3]).^2), sum((vertices_for_f[3]-vertices_for_f[1]).^2), sum((vertices_for_g[1]-vertices_for_g[2]).^2), sum((vertices_for_g[2]-vertices_for_g[3]).^2), sum((vertices_for_g[3]-vertices_for_g[1]).^2)))
                intersection=triangle_intersection(vertices_for_f, vertices_for_g)
                if(length(intersection)>=2) && (cutoffmap==nothing || cutoffmap(intersection[1]) && cutoffmap(intersection[2]))
                    push!(lines, [Point3f(intersection[1]), Point3f(intersection[2])])
                end
            end
        end
    end


    try
        foreach(line->linesegments!(ax, line; color=color, transparency=transparency, linewidth=linewidth, kwargs...), lines)
    catch e
        println("No curve in OpenGL-Mode detected! Check for the relative generality of the implicit surfaces! Error code: $(e)")
    end
end

"""
Plots a space curve, implicitly defined by two equations.
"""
function plot_implicit_curve(
    f::Function,
    g::Function;
    x_min = -2.0,
    x_max = 2.0,
    y_min = x_min,
    y_max = x_max,
    z_min = x_min,
    z_max = x_max,
    xlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (x_min, x_max),
    ylims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (y_min, y_max),
    zlims::Tuple{Union{Float64,Int},Union{Float64,Int}} = (z_min, z_max),
    resolution=(820,800),
    aspect=(1.,1.,1.),
    fontsize=17,
    show_axis=true,
    kwargs...
)
    global fig = Figure(size=resolution)
    ax = Axis3(fig[1,1], aspect=aspect)
    if !show_axis
        hidedecorations!(ax)
        hidespines!(ax)
    end
    plot_implicit_curve!(ax, f, g; xlims=xlims, ylims=ylims, zlims=zlims, kwargs...)
    return fig
end

"""
Implementation of Newton's method
"""
function newtoncorrect(equations, variables, jac, point; tol = 1e-14)
    _time = Base.time()
	q = Base.copy(point)
	global damping = 0.15
	while(norm(evaluate(equations, variables=>q)) > tol)
        if Base.time() - _time > 2
            throw(error("Error_time"))
        end
		J = evaluate.(jac, variables=>q)
		qnew = q - damping*pinv(J)*evaluate(equations, variables=>q)
		if norm(evaluate(equations, variables=>qnew)) < norm(evaluate(equations, variables=>q))
			global damping = damping*1.1
		else
			global damping = damping/2
		end
        if damping < 1e-15
            throw(error("Error_scale"))
        end
        q = qnew
		if damping > 1
			global damping = 1
		end
	end
	return q
end

"""
Compute the plane in which the `triangle` lies.
"""
function plane(triangle)
    n = cross(triangle[2] - triangle[1], triangle[3] - triangle[1])
    n /= norm(n)
    d = -dot(n, triangle[1])
    return n, d
end

signed_distance(n, d, p) = dot(n, p) + d

"""
Intersect a triangle with a plane.
"""
function triangle_plane_intersection(triangle, n, d)

    dist = [signed_distance(n,d,v) for v in triangle]

    # Snap tiny distances
    dist = [abs(x) < 1e-14 ? 0.0 : x for x in dist]

    pts = Vector{Vector{Float64}}()

    # vertices on plane
    for i in 1:3
        if dist[i] == 0.0
            push!(pts, triangle[i])
        end
    end

    edges = [(1,2),(2,3),(3,1)]

    for (i,j) in edges

        di = dist[i]
        dj = dist[j]

        if di*dj < 0

            t = di/(di-dj)

            p = triangle[i] + t*(triangle[j]-triangle[i])

            push!(pts,p)

        end

    end

    # remove duplicates
    uniquepts = Vector{Vector{Float64}}()
    for p in pts
        if !any(norm(p-q) < 1e-14 for q in uniquepts)
            push!(uniquepts,p)
        end
    end

    if length(uniquepts)==0
        return []
    elseif length(uniquepts)==1
        return uniquepts
    else
        return uniquepts[1:2]
    end
end

function overlap_segments(seg1, seg2)

    p1,p2 = seg1
    q1,q2 = seg2

    dir = p2-p1
    L = norm(dir)

    if L < 1e-14
        return []
    end

    dir /= L

    t1 = 0.0
    t2 = dot(p2-p1,dir)

    s1 = dot(q1-p1,dir)
    s2 = dot(q2-p1,dir)

    a,b = sort([t1,t2])
    c,d = sort([s1,s2])

    lo = max(a,c)
    hi = min(b,d)

    if hi < lo-1e-14
        return []
    elseif abs(hi-lo)<1e-14
        return [p1 + lo*dir]
    else
        return [p1 + lo*dir,
                p1 + hi*dir]
    end
end


function project2D(p, axis)
    if axis == 1
        return [p[2], p[3]]
    elseif axis == 2
        return [p[1], p[3]]
    else
        return [p[1], p[2]]
    end
end


orient(a,b,c) =
    (b[1]-a[1])*(c[2]-a[2]) -
    (b[2]-a[2])*(c[1]-a[1])


function line_intersection(a,b,c,d)

    r = b-a
    s = d-c

    denom = r[1]*s[2] - r[2]*s[1]

    abs(denom) < 1e-14 && error("Parallel lines")

    t = ((c[1]-a[1])*s[2] - (c[2]-a[2])*s[1]) / denom

    return a + t*r
end


function clip_edge(poly, A, B)
    isempty(poly) && return poly

    result = Vector{Vector{Float64}}()

    n = length(poly)

    for i in 1:n

        P = poly[i]
        Q = poly[mod1(i+1,n)]

        insideP = orient(A,B,P) >= -1e-14
        insideQ = orient(A,B,Q) >= -1e-14

        if insideP && insideQ

            push!(result,Q)

        elseif insideP && !insideQ

            push!(result,line_intersection(P,Q,A,B))

        elseif !insideP && insideQ

            push!(result,line_intersection(P,Q,A,B))
            push!(result,Q)

        end
    end

    return result

end


function unique_vertices(poly)

    out = Vector{Vector{Float64}}()

    for p in poly
        if isempty(out) || norm(p-out[end]) > 1e-14
            push!(out,p)
        end
    end

    if length(out)>1 && norm(out[1]-out[end])<1e-14
        pop!(out)
    end

    return out

end


function coplanar_triangle_intersection(tri1, tri2, normal)
    axis = argmax(abs.(normal))

    T1 = [project2D(v,axis) for v in tri1]
    T2 = [project2D(v,axis) for v in tri2]

    poly = copy(T1)

    # clip against each edge of triangle 2
    for i in 1:3
        A = T2[i]
        B = T2[mod(i,3)+1]
        poly = clip_edge(poly,A,B)
        isempty(poly) && return Vector{Vector{Float64}}()
    end

    poly = unique_vertices(poly)

    # Lift back into 3D.
    #
    # Since every output vertex is either
    #  - an original vertex, or
    #  - lies on an edge of one triangle,
    # we reconstruct the missing coordinate
    # using the plane equation.

    n = normal / norm(normal)
    d = -dot(n,tri1[1])

    out = Vector{Vector{Float64}}()

    for p in poly

        if axis == 1

            y,z = p
            x = -(n[2]*y + n[3]*z + d)/n[1]
            push!(out,[x,y,z])

        elseif axis == 2

            x,z = p
            y = -(n[1]*x + n[3]*z + d)/n[2]
            push!(out,[x,y,z])

        else

            x,y = p
            z = -(n[1]*x + n[2]*y + d)/n[3]
            push!(out,[x,y,z])

        end

    end

    return out

end


function triangle_intersection(triangle1,triangle2)
    n1,d1 = plane(triangle1)
    n2,d2 = plane(triangle2)

    # Parallel planes?
    if norm(cross(n1,n2)) < 1e-14

        # Different planes
        if abs(dot(n1,triangle2[1])+d1) > 1e-14
            return []
        end

        return coplanar_triangle_intersection(triangle1, triangle2, n1)
    end

    seg1 = triangle_plane_intersection(triangle1,n2,d2)
    seg2 = triangle_plane_intersection(triangle2,n1,d1)

    if isempty(seg1) || isempty(seg2)
        return []
    end

    if length(seg1)==1 && length(seg2)==1
        if norm(seg1[1]-seg2[1])<1e-14
            return seg1
        else
            return []
        end
    elseif length(seg1)==2 && length(seg2)==2
        return overlap_segments(seg1,seg2)
    else
        # point vs segment
        pt = length(seg1)==1 ? seg1[1] : seg2[1]
        seg = length(seg1)==2 ? seg1 : seg2

        a,b = seg
        ab = b-a
        t = dot(pt-a,ab)/dot(ab,ab)

        proj = a+t*ab

        if norm(proj-pt)<1e-14 &&
           -1e-14<=t<=1+1e-14
            return [pt]
        end

        return []
    end
end


function euler_step(jacobian, variables, step_size::Float64, prev_flex::Vector{Float64}, point::Union{Vector{Int},Vector{Float64}})
    J = evaluate(jacobian, variables=>point)
    flex_space = nullspace(J)
    flex_coefficients = pinv(flex_space) * prev_flex
    predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in 1:length(flex_coefficients))
    predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
    return point+step_size*predicted_inf_flex, predicted_inf_flex
end


end
