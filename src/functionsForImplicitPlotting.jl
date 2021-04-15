export plot_implicit_surface, plot_implicit_surface!, plot_implicit_curve, plot_implicit_curve!

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
    shading = true,
    wireframe=false,
    MarchingModeIsCubes=true
)
    try
        f([1,1,1])
    catch
        throw("The specified function does not take 3 arguments.")
    end
    implicit_mesh = GeometryBasics.Mesh(f,Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), MarchingModeIsCubes ? MarchingCubes() : MarchingTetrahedra())

    if wireframe
        GLMakie.wireframe!(scene, implicit_mesh, shading=shading, color=color, transparency=transparency)
    else
        vertices = GLMakie.decompose(Point{3, Float64}, implicit_mesh)
        triangles = GLMakie.decompose(TriangleFace{Int}, implicit_mesh)
        GLMakie.mesh!(scene, vertices, triangles, shading=shading, color=color, transparency=transparency)
    end
    return(scene)
end

function plot_implicit_surface(
    f;
    show_axis = true,
    resolution=(800,800),
    scale_plot=false,
    kwargs...
)
    scene = GLMakie.Scene(resolution=resolution, scale_plot=scale_plot, camera=cam3d!, show_axis=show_axis)
    plot_implicit_surface!(scene, f; kwargs...)
    return(scene)
end

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
    kwargs...
)
    try
        f([1,1,1])
        g([1,1,1])
    catch
        throw("One of the specified function does not take 3 arguments.")
    end
    f_implicit_mesh = GeometryBasics.Mesh(f, Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingCubes())
    g_implicit_mesh = GeometryBasics.Mesh(g, Rect(Vec(xlims[1], ylims[1],zlims[1]),
                            Vec(xlims[2]-xlims[1], ylims[2]-ylims[1],zlims[2]-zlims[1])), samples=samples, MarchingCubes())

    lines=[]
    for triangle_f in f_implicit_mesh
        for triangle_g in g_implicit_mesh
            if(min(sum((triangle_f[1]-triangle_g[1]).^2), sum((triangle_f[2]-triangle_g[2]).^2), sum((triangle_f[3]-triangle_g[3]).^2))
                    < max(sum((triangle_f[1]-triangle_f[2]).^2), sum((triangle_f[2]-triangle_f[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2), sum((triangle_g[1]-triangle_g[2]).^2), sum((triangle_g[2]-triangle_g[3]).^2), sum((triangle_g[3]-triangle_g[1]).^2)))
                intersection=intersectTwoTriangles(triangle_f, triangle_g)
                if(length(intersection)>=2)
                    push!(lines, [GLMakie.Point3f0(intersection[1]), GLMakie.Point3f0(intersection[2])])
                end
            end
        end
    end
    foreach(line->linesegments!(scene, line; color=color, linewidth=linewidth, kwargs...), lines)
    try
        xlims!(scene, (xlims[1],xlims[2]))
        ylims!(scene, (ylims[1],ylims[2]))
        zlims!(scene, (zlims[1],zlims[2]))
    catch e
        println("No curve was detected! Check for the relative generality of the implicit surfaces!")
    end
    return(scene)
end

function plot_implicit_curve(
    f,
    g;
    resolution=(800,800),
    scale_plot=false,
    show_axis=true,
    kwargs...
)
    scene = Scene(resolution=resolution, scale_plot=scale_plot, camera=cam3d!, show_axis=show_axis)
    plot_implicit_curve!(scene, f, g; kwargs...)
    return(scene)
end

function intersectTwoTriangles(triangle1, triangle2)
    points1 = [[point[1],point[2],point[3]] for point in triangle1]
    points2 = [[point[1],point[2],point[3]] for point in triangle2]
    polyhedron1 = polyhedron(vrep(points1))
    polyhedron2 = polyhedron(vrep(points2))
    intersection = intersect(polyhedron1, polyhedron2)
    return(vrep(intersection).points.points)
end
