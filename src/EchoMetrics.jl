module EchoMetrics

# Urmy, S. S., Horne, J. K., & Barbee, D. H. (2012). Measuring the
# vertical distributional variability of pelagic fauna in Monterey
# Bay. ICES Journal of Marine Science, 69(2), 184-196.

using Statistics
using EchogramUtils

#export density, abundance, location, dispersion, occupied_area, evenness, aggregation

export mean_volume_backscattering_strength, area_backscattering_strength, centre_of_mass,
    inertia, proportion_occupied,  equivalent_area, index_of_aggregation

"""
    mean_volume_backscattering_strength(Sv, H)

MacLennan et al. (2002).
"""
function mean_volume_backscattering_strength(Sv, H) # density
    sv = db2pow.(skipmissing(Sv))
    pow2db(sum(sv) / H)
end

"""
    area_backscattering_strength(Sv)

MacLennan et al. (2002).
"""
function area_backscattering_strength(Sv) # abundance
    sv = db2pow.(skipmissing(Sv))
    pow2db(sum(sv))
end

"""
    centre_of_mass(Sv,z)

Bez and Rivoirard (2001); Woillez et al. (2007)
"""
function centre_of_mass(Sv,z) # location
    sv = db2pow.(Sv)
    sum(z .* sv) / sum(sv)
end

"""
Bez and Rivoirard (2001); Woillez et al. (2007)
"""             
function inertia(Sv, z, CM)
    sv = db2pow.(Sv)
    sum((z .- CM).^2 .* sv) / sum(sv)
end

"""
    inertia(Sv, z)

"""
function inertia(Sv, z) 
    CM= centre_of_mass(Sv,z)
    inertia(Sv, z, CM)
end

"""
    proportion_occupied(Sv,thresh)

"""
function proportion_occupied(Sv,thresh)
    sum(Sv .> thresh) / length(Sv)
end

"""
     equivalent_area(Sv)

"""
function equivalent_area(Sv) # evenness
    sv = db2pow.(Sv)
    sum(sv)^2 / sum(sv.^2)
end

"""
    index_of_aggregation(Sv)


"""
function index_of_aggregation(Sv) # aggregation
    1 / equivalent_area(Sv)
end

end # module
