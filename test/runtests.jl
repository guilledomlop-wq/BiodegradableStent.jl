# 1. ACTIVATE THE LOCAL ENVIRONMENT
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using SparseArrays

# 2. LOAD THE PACKAGE NATIVELY
using BiodegradableStent

@testset "BiodegradableStent.jl Test Suite" begin
    
    @testset "Materials & Multiple Dispatch" begin
        # Verify material instantiation and type stability (promotion)
        # Mixed types (Int and Float) should be promoted automatically to Float64
        custom_mg = MagnesiumWE43("Test WE43", 44e9, 0.28, 1840, 250e6, 2.5e9)
        @test custom_mg.E isa Float64
        @test custom_mg.rho == 1840.0
        
        custom_zn = ZincPure("Test Zn", 90e9, 0.25, 7140.0, 120e6, 800e6)
        @test custom_zn.name == "Test Zn"
    end

    @testset "Mesh Generation & Topology" begin
        # Cross-sectional properties check
        w_width, t_depth = 0.0001, 0.0001 # 100 microns
        A, I = get_cross_section_props(w_width, t_depth)
        
        @test A == 1e-8
        @test I > 0.0

        # Auxetic Mesh Generation (2x2 cells)
        nodes, elements, w_arr, A_arr, I_arr = generate_mesh_topology(2, 2, 0.002, 0.002, -15.0, w_width, A, I)
        
        @test !isempty(nodes)
        @test !isempty(elements)
        # Verify that property arrays are correctly pre-allocated for HPC
        @test length(w_arr) == length(elements)
        @test length(A_arr) == length(elements)
        
        # Clinical metrics check
        vol, mass = calculate_stent_mass(nodes, elements, A_arr, WE43.rho)
        @test vol > 0.0
        @test mass > 0.0
    end

    @testset "Linear FEM Assembly (SparseArrays)" begin
        # 1. Local Stiffness Matrix Validation
        E_mod, A_area, I_inertia, L_length = 44e9, 1e-8, 8.3e-18, 0.002
        k_local = build_local_stiffness(E_mod, A_area, I_inertia, L_length)
        @test size(k_local) == (6, 6)
        
        # 2. Global Assembly & Sparsity Check
        nodes, elements, _, A_arr, I_arr = generate_mesh_topology(2, 2, 0.002, 0.002, -15.0, 0.0001, 1e-8, 8.3e-18)
        E_arr = fill(44e9, length(elements))
        
        K_global = assemble_global_matrix(nodes, elements, E_arr, A_arr, I_arr)
        
        # Verify it is an optimized SparseMatrixCSC
        @test K_global isa SparseMatrixCSC
        # Verify dimensions (3 DOFs per node: u, v, theta)
        @test size(K_global) == (length(nodes) * 3, length(nodes) * 3)
        # Verify it successfully populated entries
        @test nnz(K_global) > 0
    end

    @testset "Degradation: Shrinking Core Kinetics" begin
        # Verify uniform mass loss algorithms and physical boundaries
        w_initial = [0.0001, 0.0001] # 100 microns
        delta_days = 30.0
        corrosion_rate = 0.02 # mm/year
        
        # 1. Standard Degradation (Magnesium)
        w_deg_mg = calculate_corrosion_step(WE43, copy(w_initial), delta_days, corrosion_rate)
        @test w_deg_mg[1] < w_initial[1]
        
        # 2. Standard Degradation (Zinc)
        w_deg_zn = calculate_corrosion_step(PureZinc, copy(w_initial), delta_days, corrosion_rate)
        @test w_deg_zn[1] < w_initial[1]
        
        # 3. Stiffness Floor Protection (Extreme Time)
        # Ensures that a strut never reaches a negative thickness, preventing mathematical singularities
        extreme_time_days = 10000.0
        w_deg_extreme = calculate_corrosion_step(WE43, copy(w_initial), extreme_time_days, corrosion_rate)
        @test w_deg_extreme[1] == 0.0 
    end

end