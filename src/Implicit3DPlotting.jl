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

import GLMakie: xlims!, ylims!, zlims!, wireframe!, linesegments!, mesh!, Scene, cam3d!, Point3f0, scatter!, scatter, scale!, plot
import GLMakie as GLMakiePlottingLibrary
import WGLMakie as WGLMakiePlottingLibrary
import Meshing: MarchingCubes, MarchingTetrahedra
import GLMakie.GeometryBasics: Mesh, Rect, Vec, decompose, TriangleFace, Point, SimpleFaceView
import Polyhedra: vrep, intersect, polyhedron

export plot_implicit_surface,
       plot_implicit_surface!,
       plot_implicit_curve,
       plot_implicit_curve!,
#INFO: The following packages are not maintained by me. Find them here: https://github.com/MakieOrg/Makie.jl
       GLMakiePlottingLibrary,
       WGLMakiePlottingLibrary


function plot_implicit_surface!(
    fig::GLMakiePlottingLibrary.Figure,
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
    ax = fig[1,1]

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
    resolution=(800,800),
    aspect=(1.,1.,1.),
    WGLMode=false,
    in_line=false,
    transparency=true, 
    kwargs...
)
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
        global fig = WGLMakiePlottingLibrary.Scene(resolution=resolution, scale_plot=scale_plot, camera=WGLMakiePlottingLibrary.cam3d!, show_axis=show_axis)
    else
        GLMakiePlottingLibrary.activate!()
        #GLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        global fig = GLMakiePlottingLibrary.Figure(size=resolution)
        ax = GLMakiePlottingLibrary.Axis3(fig[1,1], aspect = aspect)
        if !show_axis
            GLMakiePlottingLibrary.hidespines!(ax)
            GLMakiePlottingLibrary.hidedecorations!(ax)
        end
    end
    plot_implicit_surface!(fig, f; WGLMode=WGLMode, transparency=transparency, kwargs...)
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
    samples::Tuple{Int,Int,Int}=(30,30,30),
    linewidth::Union{Float64,Int}=1.5,
    MarchingModeIsCubes::Bool = true,
    WGLMode::Bool = false,
    transparency::Bool = true,
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
        g([1,1,1])
    catch
        throw("One of the specified function does not take 3 arguments.")
    end
    f_implicit_mesh = Mesh(f, Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())
    filter!(triang->cutoffmap(triang[1]) && cutoffmap(triang[2]) && cutoffmap(triang[3]), f_implicit_mesh)
    g_implicit_mesh = Mesh(g, Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())
    filter!(triang->cutoffmap(triang[1]) && cutoffmap(triang[2]) && cutoffmap(triang[3]), g_implicit_mesh)

    lines=[]
    for triangle_f in f_implicit_mesh, triangle_g in g_implicit_mesh
        if(min(sum((triangle_f[1]-triangle_g[1]).^2), sum((triangle_f[2]-triangle_g[2]).^2), sum((triangle_f[3]-triangle_g[3]).^2))
                < max(sum((triangle_f[1]-triangle_f[2]).^2), sum((triangle_f[2]-triangle_f[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2), sum((triangle_g[1]-triangle_g[2]).^2), sum((triangle_g[2]-triangle_g[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2)))
            intersection=intersectTwoTriangles(triangle_f, triangle_g)
            if(length(intersection)>=2)
                push!(lines, [Point3f0(intersection[1]), Point3f0(intersection[2])])
            end
        end
    end

    ax = fig[1,1]

    if WGLMode
        # 3*linewidth, as the lines seem to be drawn way thinner than in GLMakie
        try
            foreach(line->WGLMakiePlottingLibrary.linesegments!(ax, line; transparency=transparency, color=color, linewidth=3*linewidth, kwargs...), lines)
            WGLMakiePlottingLibrary.xlims!(ax, (xlims[1],xlims[2]))
            WGLMakiePlottingLibrary.ylims!(ax, (ylims[1],ylims[2]))
            WGLMakiePlottingLibrary.zlims!(ax, (zlims[1],zlims[2]))
        catch e
            println("No curve in WebGL-Mode detected! Check for the relative generality of the implicit surfaces! Error code: $(e)")
        end
    else
        try
            foreach(line->linesegments!(ax, line; color=color, transparency=transparency, linewidth=linewidth, kwargs...), lines)
            xlims!(ax, (xlims[1],xlims[2]))
            ylims!(ax, (ylims[1],ylims[2]))
            zlims!(ax, (zlims[1],zlims[2]))
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
    resolution=(800,800),
    aspect=(1.,1.,1.),
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
        ax = GLMakiePlottingLibrary.Axis3(fig[1,1], aspect = aspect)
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

end
