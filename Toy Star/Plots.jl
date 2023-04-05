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