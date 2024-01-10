module Implicit3DPlotting

export plot_implicit_surface,
       plot_implicit_surface!,
       plot_implicit_curve,
       plot_implicit_curve!,
       GLMakiePlottingLibrary,
       WGLMakiePlottingLibrary

import GLMakie: xlims!, ylims!, zlims!, wireframe!, linesegments!, mesh!, Scene, cam3d!, Point3f0, scatter!, scatter, scale!
import GLMakie as GLMakiePlottingLibrary
import WGLMakie as WGLMakiePlottingLibrary
import Meshing: MarchingCubes, MarchingTetrahedra
import GLMakie.GeometryBasics: Mesh, Rect, Vec, decompose, TriangleFace, Point
import Polyhedra: vrep, intersect, polyhedron

#TODO: add shading contour function with colormap=:viridis,

"""
Adds an implicitly defined surface to the scene.
"""
function plot_implicit_surface!(
    ax,
    f;
    x_min = -3.0,
    xmin = x_min,
    x_max = 3.0,
    xmax = x_max,
    y_min = xmin,
    ymin = y_min,
    y_max = xmax,
    ymax = y_max,
    z_min = xmin,
    zmin = z_min,
    z_max = xmax,
    zmax = z_max,
    xlims = (xmin, xmax),
    ylims = (ymin, ymax),
    zlims = (zmin, zmax),
    color = :steelblue,
    transparency = true,
    samples=(35,35,35),
    wireframe=false,
    MarchingModeIsCubes=true,
    WGLMode = false,
    zcolormap=nothing,
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
    if wireframe
        if WGLMode
            WGLMakiePlottingLibrary.wireframe!(ax, implicit_mesh, color=color, transparency=transparency)
        else
            wireframe!(ax, implicit_mesh, color=color, transparency=transparency)
        end
    else
        vertices = decompose(Point{3, Float64}, implicit_mesh)
        triangles = decompose(TriangleFace{Int}, implicit_mesh)
        if WGLMode
            if isnothing(zcolormap)
                WGLMakiePlottingLibrary.mesh!(ax, vertices, triangles, color=color, transparency=transparency, kwargs...)
            else
                colors = [v[3] for v in vertices]
                WGLMakiePlottingLibrary.mesh!(ax, vertices, triangles, color=colors, transparency=transparency, colormap=zcolormap, kwargs...)
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

    return(ax)
end

"""
Plots an implicitly defined surface.
"""
function plot_implicit_surface(
    f;
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
        #WGLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        scene = WGLMakiePlottingLibrary.Scene(resolution=resolution, scale_plot=scale_plot, camera=WGLMakiePlottingLibrary.cam3d!, show_axis=show_axis)
    else
        GLMakiePlottingLibrary.activate!()
        #GLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        scene = GLMakiePlottingLibrary.Figure(size=resolution)
        ax = GLMakiePlottingLibrary.Axis3(scene[1,1], aspect = aspect)
        if !show_axis
            GLMakiePlottingLibrary.hidespines!(ax)
            GLMakiePlottingLibrary.hidedecorations!(ax)
        end
    end
    plot_implicit_surface!(ax, f; WGLMode=WGLMode, transparency=transparency, kwargs...)
    return(scene)
end

"""
Adds a space curve, implicitly defined by two equations, to a scene.
"""
function plot_implicit_curve!(
    ax,
    f,
    g;
    x_min = -2.0,
    xmin = x_min,
    x_max = 2.0,
    xmax = x_max,
    y_min = xmin,
    ymin = y_min,
    y_max = xmax,
    ymax = y_max,
    z_min = xmin,
    zmin = z_min,
    z_max = xmax,
    zmax = z_max,
    xlims = (xmin, xmax),
    ylims = (ymin, ymax),
    zlims = (zmin, zmax),
    color = :steelblue,
    samples=(30,30,30),
    linewidth=1.5,
    MarchingModeIsCubes=true,
    WGLMode = false,
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
    g_implicit_mesh = Mesh(g, Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())

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
    if WGLMode
        # 3*linewidth, as the lines seem to be drawn way thinner than in GLMakie
        try
            foreach(line->WGLMakiePlottingLibrary.linesegments!(ax, line; color=color, linewidth=3*linewidth, kwargs...), lines)
            WGLMakiePlottingLibrary.xlims!(ax, (xlims[1],xlims[2]))
            WGLMakiePlottingLibrary.ylims!(ax, (ylims[1],ylims[2]))
            WGLMakiePlottingLibrary.zlims!(ax, (zlims[1],zlims[2]))
        catch e
            println("No curve in WebGL-Mode detected! Check for the relative generality of the implicit surfaces!")
        end
    else
        try
            foreach(line->linesegments!(ax, line; color=color, linewidth=linewidth, kwargs...), lines)
            GLMakiePlottingLibrary.xlims!(ax, (xlims[1],xlims[2]))
            GLMakiePlottingLibrary.ylims!(ax, (ylims[1],ylims[2]))
            GLMakiePlottingLibrary.zlims!(ax, (zlims[1],zlims[2]))
        catch e
            println("No curve in OpenGL-Mode detected! Check for the relative generality of the implicit surfaces!")
        end
    end

    return(ax)
end

"""
Plots a space curve, implicitly defined by two equations.
"""
function plot_implicit_curve(
    f,
    g;
    resolution=(800,800),
    aspect=(1.,1.,1.),
    show_axis=true,
    in_line=false,
    WGLMode=false,
    kwargs...
)
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
        #WGLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        scene = WGLMakiePlottingLibrary.Scene(resolution=resolution, camera=WGLMakiePlottingLibrary.cam3d!, show_axis=show_axis)
    else
        GLMakiePlottingLibrary.activate!()
        #GLMakiePlottingLibrary.AbstractPlotting.inline!(in_line)
        scene = GLMakiePlottingLibrary.Figure(size=resolution)
        ax = GLMakiePlottingLibrary.Axis3(figure[1,1], aspect = aspect)
        if !show_axis
            GLMakiePlottingLibrary.hidespines!(ax)
            GLMakiePlottingLibrary.hidedecorations!(ax)
        end
    end
    plot_implicit_curve!(ax, f, g; WGLMode=WGLMode, kwargs...)

    return(scene)
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
