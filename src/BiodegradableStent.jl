module BiodegradableStent

# Include internal submodules
include("Materials.jl")
include("MeshGenerator.jl")
include("LinearFEM.jl")
include("Degradation.jl")

# Expose submodules internally for re-exporting the API
using .Materials
using .MeshGenerator
using .LinearFEM
using .Degradation

# =========================================================================
# PUBLIC API EXPORTS
# =========================================================================

# From Materials.jl
export AbstractBioMaterial, MagnesiumWE43, ZincPure
export WE43, PureZinc

# From MeshGenerator.jl
export generate_mesh_topology, get_cross_section_props, calculate_stent_mass

# From LinearFEM.jl
export build_local_stiffness, assemble_global_matrix, solve_elastic_system

# From Degradation.jl
export calculate_corrosion_step

end