module Materials

export AbstractBioMaterial, MagnesiumWE43, ZincPure
export WE43, PureZinc

# 0. ABSTRACT TYPE HIERARCHY
"""
Abstract base type for all biodegradable stent materials.
Enables multiple dispatch for degradation and mechanical behavior.
"""
abstract type AbstractBioMaterial{T<:Real} end

# 1. MATERIAL PROPERTIES
"""
MagnesiumWE43{T<:Real}

Defines the mechanical and physical properties of a bioabsorbable material.
Parametrized with `T<:Real` to support ForwardDiff (dual numbers) for future optimization.
"""
struct MagnesiumWE43{T} <: AbstractBioMaterial{T}
    name::String
    E::T                 # Young's Modulus (Pa)
    nu::T                # Poisson ratio
    rho::T               # Density (kg/m^3)
    yield_strength::T    # Yield strength (Pa)
    hardening_modulus::T # Tangent hardening modulus (Pa)

    # Inner Constructor: Enforces strict type stability for peak CPU performance.
    # Ensures all numerical fields share the exact same type 'T' in memory.
    function MagnesiumWE43(name::String, E::T, nu::T, rho::T, 
                           ys::T, hm::T) where T<:Real
        new{T}(name, E, nu, rho, ys, hm)
    end
end

# Outer Constructor: Facilitates instantiation by automatically promoting mixed 
# input types (e.g., Int and Float64) to a common, optimal type 'T'.
function MagnesiumWE43(name::String, E, nu, rho, ys, hm)
    T = promote_type(typeof(E), typeof(nu), typeof(rho), typeof(ys), typeof(hm))
    return MagnesiumWE43(name, T(E), T(nu), T(rho), T(ys), T(hm))
end

# 1.1 MATERIAL PROPERTIES (ZINC)
struct ZincPure{T} <: AbstractBioMaterial{T}
    name::String
    E::T
    nu::T
    rho::T
    yield_strength::T
    hardening_modulus::T
    
    function ZincPure(name::String, E::T, nu::T, rho::T, ys::T, hm::T) where T<:Real
        new{T}(name, E, nu, rho, ys, hm)
    end
end

function ZincPure(name::String, E, nu, rho, ys, hm)
    T = promote_type(typeof(E), typeof(nu), typeof(rho), typeof(ys), typeof(hm))
    return ZincPure(name, T(E), T(nu), T(rho), T(ys), T(hm))
end

# 2. REFERENCE MATERIALS
# Declared as 'const' to prevent type changes and optimize memory allocation.

const WE43 = MagnesiumWE43(
    "Magnesium WE43",
    44e9,      # E = 44 GPa
    0.28,      # nu
    1840.0,    # rho
    250e6,     # yield_strength (250 MPa)
    2.5e9      # hardening_modulus
)

const PureZinc = ZincPure(
    "Pure Zinc (Zn)",
    90e9,      # E = 90 GPa
    0.25,      # nu
    7140.0,    # rho
    120e6,     # yield_strength (120 MPa)
    800e6      # hardening_modulus 
)

# Note on Degradation Protocol:
# To maintain Data-Oriented Design and stack allocation efficiency, this struct remains immutable.
# Material degradation will be simulated by creating new instances of BioMaterial rather than mutating.

end