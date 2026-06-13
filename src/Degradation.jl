module Degradation

export calculate_corrosion_step

using ..Materials: AbstractBioMaterial, MagnesiumWE43, ZincPure

# =========================================================================
# THE MAIN INTERFACE (Public Interface)
# Uses Multiple Dispatch to route arrays to the correct chemical recipe
# =========================================================================
"""
calculate_corrosion_step(material, w_array, delta_t, P_rate)

Public Interface for basic chemical degradation.
Routes the geometry arrays to the correct kinetic model based on the Material type.
"""

# Dispatch 1: Magnesium (Base uniform degradation)
function calculate_corrosion_step(material::MagnesiumWE43{T}, w_array::Vector{T}, delta_t::T, P_rate::T) where T<:Real
    return degrade_magnesium_uniform(material, w_array, delta_t, P_rate)
end

# Dispatch 2: Zinc (Alternative degradation profile)
function calculate_corrosion_step(material::ZincPure{T}, w_array::Vector{T}, delta_t::T, P_rate::T) where T<:Real
    return degrade_zinc_uniform(material, w_array, delta_t, P_rate)
end

# =========================================================================
# SPECIFIC RECIPES (Private Kinetic Models)
# =========================================================================
"""
degrade_magnesium_uniform(material, w_array, delta_t, P_rate)

Calculates the residual thickness of a strut array based on the uniform 
corrosion model proposed by Yufeng Zheng (Shrinking Core Model).

Inputs:
- w_array: Current thickness/width vector (m)
- delta_t: Time step duration (days)
- P_rate: Base corrosion rate in mm/year (Standard clinical metric for WE43)
"""
function degrade_magnesium_uniform(material::MagnesiumWE43{T}, w_array::Vector{T}, delta_t::T, P_rate::T) where T<:Real
    num_elements = length(w_array)
    w_new = Vector{T}(undef, num_elements)
    
    # 1. Unit conversion: mm/year -> meters/day
    # (P_rate / 1000.0) converts to m/year
    # Dividing by 365.25 converts to m/day
    v_corr_base = (P_rate / T(1000.0)) / T(365.25)
    
    for i in 1:num_elements
        # 2. Shrinking Core Model application:
        # Chemical attack occurs uniformly on both exposed faces of the strut.
        # Therefore, the width reduction is twice the corrosion velocity.
        # Equation: w(t) = w(t-1) - 2 * v_corr * delta_t
        w_degraded = w_array[i] - (T(2.0) * v_corr_base * delta_t)
        
        # 3. Stiffness floor: Material thickness cannot physically be negative.
        # We maintain the strict floor threshold here.
        w_new[i] = max(w_degraded, T(0.0)) 
    end
    
    return w_new
end

"""
degrade_zinc_uniform(material, w_array, delta_t, P_rate)

Zinc degradation involves distinct oxide layer formation, approximated here
via the uniform shrinking core approach for baseline comparisons.
"""
function degrade_zinc_uniform(material::ZincPure{T}, w_array::Vector{T}, delta_t::T, P_rate::T) where T<:Real
    num_elements = length(w_array)
    w_new = Vector{T}(undef, num_elements)
    v_corr_base = (P_rate / T(1000.0)) / T(365.25)
    
    for i in 1:num_elements
        # Pure uniform degradation for the baseline portfolio model
        w_degraded = w_array[i] - (T(2.0) * v_corr_base * delta_t)
        w_new[i] = max(w_degraded, T(0.0)) 
    end
    
    return w_new
end

end