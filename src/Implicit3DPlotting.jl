module Implicit3DPlotting

export plot_implicit_surface,
       plot_implicit_surface!,
       plot_implicit_curve,
       plot_implicit_curve!,
       GLMakiePlottingLibrary,
       WGLMakiePlottingLibrary

import GLMakie: xlims!, ylims!, zlims!, wireframe!, linesegments!, mesh!, Scene, cam3d!, Point3f0, scatter!, scatter
import GLMakie as GLMakiePlottingLibrary
import WGLMakie as WGLMakiePlottingLibrary
import Meshing: MarchingCubes, MarchingTetrahedra
import GLMakie.GeometryBasics: Mesh, Rect, Vec, decompose, TriangleFace, Point
import Polyhedra: vrep, intersect, polyhedron

"""
Adds an implicitly defined surface to the scene.
"""
function plot_implicit_surface!(
    scene,
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
    shading = true,
    wireframe=false,
    MarchingModeIsCubes=true,
    WGLMode = false,
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
            WGLMakiePlottingLibrary.wireframe!(scene, implicit_mesh, shading=shading, color=color, transparency=transparency)
        else
            wireframe!(scene, implicit_mesh, shading=shading, color=color, transparency=transparency)
        end
    else
        vertices = decompose(Point{3, Float64}, implicit_mesh)
        triangles = decompose(TriangleFace{Int}, implicit_mesh)
        if WGLMode
            WGLMakiePlottingLibrary.mesh!(scene, vertices, triangles, shading=shading, color=color, transparency=transparency)
        else
            mesh!(scene, vertices, triangles, shading=shading, color=color, transparency=transparency)
        end
    end

    return(scene)
end

"""
Plots an implicitly defined surface.
"""
function plot_implicit_surface(
    f;
    show_axis = true,
    resolution=(800,800),
    scale_plot=false,
    WGLMode=false,
    in_line=false,
    kwargs...
)
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
        scene = WGLMakiePlottingLibrary.Scene(resolution=resolution, scale_plot=scale_plot, camera=WGLMakiePlottingLibrary.cam3d!, show_axis=show_axis)
    else
        GLMakiePlottingLibrary.activate!()
        scene = Scene(resolution=resolution, scale_plot=scale_plot, camera=cam3d!, show_axis=show_axis)
    end
    plot_implicit_surface!(scene, f; WGLMode=WGLMode, kwargs...)
    return(scene)
end

"""
Adds a space curve, implicitly defined by two equations, to a scene.
"""
function plot_implicit_curve!(
    scene,
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
    for triangle_f in f_implicit_mesh
        for triangle_g in g_implicit_mesh
            if(min(sum((triangle_f[1]-triangle_g[1]).^2), sum((triangle_f[2]-triangle_g[2]).^2), sum((triangle_f[3]-triangle_g[3]).^2))
                    < max(sum((triangle_f[1]-triangle_f[2]).^2), sum((triangle_f[2]-triangle_f[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2), sum((triangle_g[1]-triangle_g[2]).^2), sum((triangle_g[2]-triangle_g[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2)))
                intersection=intersectTwoTriangles(triangle_f, triangle_g)
                if(length(intersection)>=2)
                    push!(lines, [Point3f0(intersection[1]), Point3f0(intersection[2])])
                end
            end
        end
    end
    if WGLMode
        # 3*linewidth, as the lines seem to be drawn way thinner than in GLMakie
        foreach(line->WGLMakiePlottingLibrary.linesegments!(scene, line; color=color, linewidth=3*linewidth, kwargs...), lines)
        try
            WGLMakiePlottingLibrary.xlims!(scene, (xlims[1],xlims[2]))
            WGLMakiePlottingLibrary.ylims!(scene, (ylims[1],ylims[2]))
            WGLMakiePlottingLibrary.zlims!(scene, (zlims[1],zlims[2]))
        catch e
            println("No curve in WebGL-Mode detected! Check for the relative generality of the implicit surfaces!")
        end
    else
        foreach(line->linesegments!(scene, line; color=color, linewidth=linewidth, kwargs...), lines)
        try
            xlims!(scene, (xlims[1],xlims[2]))
            ylims!(scene, (ylims[1],ylims[2]))
            zlims!(scene, (zlims[1],zlims[2]))
        catch e
            println("No curve in OpenGL-Mode detected! Check for the relative generality of the implicit surfaces!")
        end
    end
    return(scene)
end

"""
Plots a space curve, implicitly defined by two equations.
"""
function plot_implicit_curve(
    f,
    g;
    resolution=(800,800),
    scale_plot=false,
    show_axis=true,
    WGLMode=false,
    kwargs...
)
    if WGLMode
        WGLMakiePlottingLibrary.activate!()
        scene = WGLMakiePlottingLibrary.Scene(resolution=resolution, scale_plot=scale_plot, camera=WGLMakiePlottingLibrary.cam3d!, show_axis=show_axis)
    else
        GLMakiePlottingLibrary.activate!()
        scene = Scene(resolution=resolution, scale_plot=scale_plot, camera=cam3d!, show_axis=show_axis)
    end
    plot_implicit_curve!(scene, f, g; WGLMode=WGLMode, kwargs...)
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
