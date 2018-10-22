module EchoMetrics

# Urmy, S. S., Horne, J. K., & Barbee, D. H. (2012). Measuring the
# vertical distributional variability of pelagic fauna in Monterey
# Bay. ICES Journal of Marine Science, 69(2), 184-196.

using Statistics
using EchogramUtils

export density, abundance, location, dispersion, occupied_area, evenness, aggregation

"""
MacLennan et al. (2002).
"""
function mean_volume_backscattering_strength(sv, H) # density
    pow2db(sum(sv) / H)
end

"""
MacLennan et al. (2002).
"""
function area_backscattering_strength(sv) # abundance
    pow2db(sum(sv))
end

"""
Bez and Rivoirard (2001); Woillez et al. (2007)
"""
function centre_of_mass(z, sv) # location
    sum(z .* sv) / sum(sv)
end

"""
Bez and Rivoirard (2001); Woillez et al. (2007)
"""             
function inertia(z, sv, CM)
    @info("inertia")
    sum((z .- CM).^2 .* sv) / sum(sv)
end

function inertia(z, sv) 
    CM= centre_of_mass(z, sv)
    inertia(z,sv, CM)
end

function proportion_occupied(sv, z, thresh)
    sum(sv .> thresh) / length(z)
end

function equivalent_area(sv) # evenness
    sum(sv)^2 / sum(sv.^2)
end


function index_of_aggregation(sv) # aggregation
    1 / equivalent_area(sv)
end


function number_of_layers(x) # layer_structure
end

function cols(A)
    m,n = size(A)
    [A[:,i] for i in 1:n]
end

function density(A, H)
    sv = db2pow.(A)
    svs = cols(sv)

    mean_volume_backscattering_strength.(svs, H)
end

function abundance(A)
    sv = db2pow.(A)
    svs = cols(sv)
    area_backscattering_strength.(svs)
end

function location(A, R)
    sv = db2pow.(A)
    svs = cols(sv)
    rs = cols(R)
    centre_of_mass.(rs, svs)
end

function dispersion(A, R)
    sv = db2pow.(A)
    svs = cols(sv)
    zs = cols(R)
    CMs = centre_of_mass.(zs, svs)

    inertia.(zs, svs, CMs)
end

function occupied_area(A, R)
    sv = db2pow.(A)
    svs = cols(sv)
    zs = cols(R)
    thresh = db2pow(-90)
    proportion_occupied.(svs, zs, thresh)
end

function evenness(A)
    sv = db2pow.(A)
    svs = cols(sv)
    equivalent_area.(svs)
end

function aggregation(A)
    sv = db2pow.(A)
    svs = cols(sv)
    index_of_aggregation.(svs)
end

end # module
