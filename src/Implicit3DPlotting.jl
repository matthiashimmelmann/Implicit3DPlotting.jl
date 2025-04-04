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

import GLMakie: xlims!, ylims!, zlims!, wireframe!, linesegments!, mesh!, Scene, cam3d!, Point3f, scatter!, scatter, scale!, plot, DirectionalLight, LScene, RGBf, Vec3f
import GLMakie as GLMakiePlottingLibrary
import WGLMakie as WGLMakiePlottingLibrary
import Meshing: MarchingCubes, MarchingTetrahedra
import GLMakie.GeometryBasics: Mesh, Rect, Vec, decompose, TriangleFace, Point, SimpleFaceView
import Polyhedra: vrep, intersect, polyhedron
import LinearAlgebra: pinv, norm, nullspace
import HomotopyContinuation: evaluate, differentiate, @var, Expression
import IterTools: product

export plot_implicit_surface,
       plot_implicit_surface!,
       plot_implicit_curve,
       plot_implicit_curve!,
#INFO: The following packages are not maintained by me. Find them here: https://github.com/MakieOrg/Makie.jl
       GLMakiePlottingLibrary,
       WGLMakiePlottingLibrary


function plot_implicit_surface!(
    ax,
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
    color::Symbol = :steelblue,
    transparency::Bool = true,
    samples::Tuple{Int,Int,Int}=(35,35,35),
    wireframe::Bool=false,
    MarchingModeIsCubes::Bool=true,
    WGLMode::Bool = false,
    zcolormap=nothing,
    cutoffmap=nothing,
    kwargs...
    )
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
    else
        GLMakiePlottingLibrary.activate!()
    end

    try
        f([1,1,1])
    catch
        throw("The specified function does not take 3 arguments.")
    end
    implicit_mesh = Mesh(f,Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())
    
    vertices = decompose(Point{3, Float64}, implicit_mesh)
    triangles = decompose(TriangleFace{Int}, implicit_mesh)
    i = 1
    if cutoffmap!=nothing
        while i <= length(triangles)
            if cutoffmap(vertices[triangles[i][1]]) && cutoffmap(vertices[triangles[i][2]]) && cutoffmap(vertices[triangles[i][3]])
                deleteat!(triangles, i)
            else
                i=i+1
            end
        end
    end

    if wireframe
        implicit_mesh = Mesh(vertices, triangles)
        if WGLMode
            wireframe!(ax, implicit_mesh, color=color, transparency=transparency)
        else
            wireframe!(ax, implicit_mesh, color=color, transparency=transparency)
        end
    else
        if WGLMode
            if isnothing(zcolormap)
                mesh!(ax, vertices, triangles, color=color, transparency=transparency, kwargs...)
            else
                colors = [v[3] for v in vertices]
                mesh!(ax, vertices, triangles, color=colors, transparency=transparency, colormap=zcolormap, kwargs...)
            end
        else
            if isnothing(zcolormap)
                mesh!(ax, vertices, triangles, color=color, transparency=transparency, kwargs...)
            else
                colors = [v[3] for v in vertices]
                mesh!(ax, vertices, triangles, color=colors, transparency=transparency, colormap=zcolormap, kwargs...)
            end        
        end
    end
end

#=
Plots an implicitly defined surface.
=#
function plot_implicit_surface(
    f::Function;
    show_axis = true,
    resolution=(820,800),
    aspect=(1.,1.,1.),
    WGLMode=false,
    in_line=false,
    transparency=true,
    fontsize=17,
    lightingdirections = [Vec3f(t) for t in vcat(product(-1:2:1, -1:2:1, -1:2:1)...)],
    lightintensity=0.425,
    kwargs...
)
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
        global fig = WGLMakiePlottingLibrary.Scene(resolution=resolution, scale_plot=scale_plot, camera=WGLMakiePlottingLibrary.cam3d!, show_axis=show_axis)
        lights = [DirectionalLight(RGBf(lightintensity, lightintensity, lightintensity), vector) for vector in lightingdirections]
        ax = LScene(fig[1, 1], scenekw = (lights = lights, aspect = aspect, xlabelsize=fontsize, ylabelsize=fontsize, zlabelsize=fontsize), show_axis=show_axis)
    else
        GLMakiePlottingLibrary.activate!()
        #GLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        global fig = GLMakiePlottingLibrary.Figure(size=resolution)
        lights = [DirectionalLight(RGBf(lightintensity, lightintensity, lightintensity), vector) for vector in lightingdirections]
        ax = LScene(fig[1, 1], scenekw = (lights = lights, aspect = aspect, xlabelsize=fontsize, ylabelsize=fontsize, zlabelsize=fontsize), show_axis=show_axis)

    end
    plot_implicit_surface!(ax, f; WGLMode=WGLMode, transparency=transparency, kwargs...)
    return(fig)
end

"""
Adds a space curve, implicitly defined by two equations, to a scene.
"""
function plot_implicit_curve!(
    fig::GLMakiePlottingLibrary.Figure,
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
    color::Symbol = :steelblue,
    samples::Tuple{Int,Int,Int}=(8,8,8),
    linewidth::Union{Float64,Int}=2.25,
    MarchingModeIsCubes::Bool = true,
    WGLMode::Bool = false,
    step_size=0.025,
    transparency::Bool = true,
    cutoffmap=nothing,
    miniter=5,
    kwargs...
)
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
    else
        GLMakiePlottingLibrary.activate!()
    end

    try
        f([1,1,1])
        g([1,1,1])
    catch
        throw("One of the specified function does not take 3 arguments.")
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
        end

    catch e
        @warn "The given functions f and g could not be cast to polynomials. Error code: $(e)."
        f_implicit_mesh = Mesh(f, Rect(Vec(xlims[1], ylims[1],zlims[1]),
                                Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())
        g_implicit_mesh = Mesh(g, Rect(Vec(xlims[1], ylims[1],zlims[1]),
                                Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())

        for triangle_f in f_implicit_mesh, triangle_g in g_implicit_mesh
            if(min(sum((triangle_f[1]-triangle_g[1]).^2), sum((triangle_f[2]-triangle_g[2]).^2), sum((triangle_f[3]-triangle_g[3]).^2))
                    < max(sum((triangle_f[1]-triangle_f[2]).^2), sum((triangle_f[2]-triangle_f[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2), sum((triangle_g[1]-triangle_g[2]).^2), sum((triangle_g[2]-triangle_g[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2)))
                intersection=intersectTwoTriangles(triangle_f, triangle_g)
                if(length(intersection)>=2) && (cutoffmap==nothing || cutoffmap(intersection[1]) && cutoffmap(intersection[2]))
                    push!(lines, [Point3f(intersection[1]), Point3f(intersection[2])])
                end
            end
        end
    end

    ax = fig[1,1]

    if WGLMode
        # 3*linewidth, as the lines seem to be drawn way thinner than in GLMakie
        try
            foreach(line->WGLMakiePlottingLibrary.linesegments!(ax, line; transparency=transparency, color=color, linewidth=3*linewidth, kwargs...), lines)
            WGLMakiePlottingLibrary.xlims!(ax, xlims[1],xlims[2])
            WGLMakiePlottingLibrary.ylims!(ax, ylims[1],ylims[2])
            WGLMakiePlottingLibrary.zlims!(ax, zlims[1],zlims[2])
        catch e
            println("No curve in WebGL-Mode detected! Check for the relative generality of the implicit surfaces! Error code: $(e)")
        end
    else
        try
            foreach(line->linesegments!(ax, line; color=color, transparency=transparency, linewidth=linewidth, kwargs...), lines)
            xlims!(xlims[1],xlims[2])
            ylims!(ylims[1],ylims[2])
            zlims!(zlims[1],zlims[2])
        catch e
            println("No curve in OpenGL-Mode detected! Check for the relative generality of the implicit surfaces! Error code: $(e)")
        end
    end
end

"""
Plots a space curve, implicitly defined by two equations.
"""
function plot_implicit_curve(
    f::Function,
    g::Function;
    resolution=(820,800),
    aspect=(1.,1.,1.),
    fontsize=17,
    show_axis=true,
    in_line=false,
    WGLMode=false,
    kwargs...
)
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
        WGLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        global fig = WGLMakiePlottingLibrary.Scene(resolution=resolution, camera=WGLMakiePlottingLibrary.cam3d!, show_axis=show_axis)
    else
        GLMakiePlottingLibrary.activate!()
        #GLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        global fig = GLMakiePlottingLibrary.Figure(size=resolution)
    
        ax = GLMakiePlottingLibrary.Axis3(fig[1,1], aspect = aspect, xlabelsize=fontsize, ylabelsize=fontsize, zlabelsize=fontsize)
        if !show_axis
            GLMakiePlottingLibrary.hidespines!(ax)
            GLMakiePlottingLibrary.hidedecorations!(ax)
        end
    end
    plot_implicit_curve!(fig, f, g; WGLMode=WGLMode, kwargs...)
    return(fig)
end

"""
Calculates the intersection of two TriangleFaces.
"""
function intersectTwoTriangles(triangle1, triangle2)
    points1 = [[point[1],point[2],point[3]] for point in triangle1]
    points2 = [[point[1],point[2],point[3]] for point in triangle2]
    polyhedron1 = polyhedron(vrep(points1))
    polyhedron2 = polyhedron(vrep(points2))
    intersection = intersect(polyhedron1, polyhedron2)
    return(vrep(intersection).points.points)
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


function euler_step(jacobian, variables, step_size::Float64, prev_flex::Vector{Float64}, point::Union{Vector{Int},Vector{Float64}})
    J = evaluate(jacobian, variables=>point)
    flex_space = nullspace(J)[:,1]
    flex_coefficients = pinv(flex_space) * prev_flex
    predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in 1:length(flex_coefficients))
    predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
    return point+step_size*predicted_inf_flex, predicted_inf_flex
end


end
