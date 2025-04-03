import Implicit3DPlotting: plot_implicit_surface, plot_implicit_curve
using Test


@testset "hyperboloid" begin
    f(x) = sum(x.^2) - 1
    plot_implicit_surface(f; transparency=false, xlims=(-3,3), ylims=(-3,3), zlims=(-3,3))
end

@testset "gordianknotcurves" begin
    g1(x) = sum(x.^2) - 1
    g2(x) = 2*prod(x)+1-sum(x.^2)
    plot_implicit_curve(g1, g2; samples = (5,5,5), step_size=0.025)
end

@testset "cayleyspectahedron" begin
    h = x -> 2*prod(x) + 1 - sum(x.^2)
    plot_implicit_surface(h; cutoffmap=x->x[1]^2+x[2]^2+x[3]^2-3>=0, samples=(500,500,500),
           transparency=false, zcolormap=:viridis, xlims=(-1.25,1.25), ylims=(-1.25,1.25), zlims=(-1.25,1.25))
end