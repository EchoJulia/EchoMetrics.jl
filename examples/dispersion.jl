#!/usr/bin/env julia

using EchoMetrics
using SimradEK60
using SimradEK60TestData
using EchogramPyPlot
using PyPlot
using LaTeXStrings
using EchogramUtils


function main()
    @info("Start")
    
    filename = EK60_SAMPLE
    ps = SimradEK60.load(filename)
    ps38 = [p for p in ps if p.frequency == 38000]
    Sv38 = Sv(ps38) # Volume backscatter
    #al38 = alongshipangle(ps38) # Split beam angle
    #at38 = athwartshipangle(ps38)
    _R = R(ps38) # Range / depth
    _t = filetime(ps38) # timestamps

    Sv38s = pow2db.(vertically_smooth(db2pow.(Sv38),_R,thickness=1.43));
    mask = INmask(Sv38s, delta=10);
    
    Sv38[_R .< 10] .= -999
    Sv38[mask] .= -999
    
    d = dispersion(Sv38,_R)

    fig, axes = subplots(2,1, figsize=(6.5,8.5))

    ax = subplot(2,1,1)
    echogram(Sv38, vmin = -95, vmax=-50, range = maximum(_R), cmap=EK500)
    cb = plt[:colorbar]()
    cb[:set_label](L"$S_v$ dB re 1 m$^{-1}$")
    ylabel("Range / m")

    ax = subplot(2,1,2)
    plot(d)
    ylabel(L"I (m$^{-2}$)")
    xlabel("Ping")
    fig[:tight_layout]()

    filename = "dispersion.png"
    @info("Saving $filename ...")

    plt[:savefig](filename, dpi=600)

end


main()
