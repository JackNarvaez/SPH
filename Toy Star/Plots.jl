using Plots
using LaTeXStrings

function animation(T, record, dt, pos, Prop, Property, pmin, pmax, R, Name)
    anim = @animate for n in 1:record:T
        plot(pos[n, :, 1], pos[n, :, 2], seriestype=:scatter, aspect_ratio=:equal, cmap=:thermal, markerstrokewidth=0,
            marker_z = Prop[n,:], ms=1.5, ma=0.5, xlim=(-4/3*R, 4/3*R), ylim=(-4/3*R, 4/3*R), legend=false, grid=false, 
            colorbar=true, colorbar_title="\n"*Property, clims=(pmin,pmax), dpi=300)
        xlabel!(L"x")
        ylabel!(L"y")
        annotate!(0, 4/3*R, text("t = $(round(n*dt;digits=1))",8))
    end
    gif(anim, Name*".mp4", fps = 10)
end;

function animation2(T, record, dt, pos, Prop, Property, pmin, pmax, R, Name, N1)
    anim = @animate for n in 1:record:T
        plot(pos[n, 1:N1, 1], pos[n, 1:N1, 2], seriestype=:scatter, aspect_ratio=:equal, mc=:red, markerstrokewidth=0,
            ms=1.5, ma=0.5, xlim=(-4/3*R, 4/3*R), ylim=(-4/3*R, 4/3*R), legend=false, grid=false, 
            dpi=300)
        plot!(pos[n, N1+1:end, 1], pos[n, N1+1:end, 2], seriestype=:scatter, aspect_ratio=:equal, mc=:blue, markerstrokewidth=0,
            ms=1.5, ma=0.5, xlim=(-4/3*R, 4/3*R), ylim=(-4/3*R, 4/3*R), legend=false, grid=false, 
            dpi=300)
        xlabel!(L"x")
        ylabel!(L"y")
        annotate!(0, 4/3*R, text("t = $(round(n*dt;digits=1))",8))
    end
    gif(anim, Name*".mp4", fps = 10)
end;

function PlotDensity(radius, dens, dens_Theo, R, PropMax)
    P_r = sortperm(radius)
    plot(radius[P_r], dens[P_r],  seriestype=:scatter, ms=1.2, markerstrokewidth=0, color=:crimson, xlim=(0, R), ylim=(0, PropMax), labels=L"\rho_{SPH}(r)", grid=false, dpi=300)
    plot!(radius[P_r], dens_Theo[P_r], color=:blue, labels=L"\rho_{theo}(r)", grid=false, legend=true, dpi=300)
    xlabel!(L"r")
    ylabel!(L"\rho")
    savefig("Density.png")
end;