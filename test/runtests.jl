using Test
using EchoMetrics


Sv = [ -70, -60, - 50]
H = 30

density = mean_volume_backscattering_strength(Sv)

@test trunc(Int, density) == -54

abundance = area_backscattering_strength(Sv, H)

@test trunc(Int, abundance) == -39
