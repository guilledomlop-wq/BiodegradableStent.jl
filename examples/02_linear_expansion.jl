"""
Example 02: Linear FEM Expansion

This script demonstrates the high-performance SparseArrays assembly and 
the strictly linear elastic solver to simulate the axial expansion 
of the biodegradable stent. Visual results are saved locally.
"""

# 1. ACTIVATE AND INSTANTIATE THE LOCAL ENVIRONMENT
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# 2. LOAD THE PACKAGE NATIVELY
using BiodegradableStent
using Plots

# Force Plots (GR backend) into headless mode to prevent any window pop-ups in VS Code/REPL
ENV["GKSwstype"] = "100"

function run_linear_expansion()
    println("Starting Linear Expansion Simulation...")

    # 3. BASE MESH GENERATION (Re-using parameters from Example 1)
    rows = 4
    cols = 8
    h = 0.002
    l = 0.002
    theta_deg = -15.0
    width = 0.0001
    depth = 0.0001
    
    A, I = get_cross_section_props(width, depth)
    nodes, elements, w_array, A_array, I_array = generate_mesh_topology(
        rows, cols, h, l, theta_deg, width, A, I
    )
    
    # 4. MATERIAL SETUP (Using Magnesium WE43)
    # Fill an array with the Young's Modulus for each strut
    E_array = fill(WE43.E, length(elements))
    
    # 5. GLOBAL MATRIX ASSEMBLY (HPC Showcase)
    println("Assembling global stiffness matrix (SparseArrays)...")
    K_global = assemble_global_matrix(nodes, elements, E_array, A_array, I_array)
    println("Matrix assembled successfully. Total DOFs: ", size(K_global, 1))

    # 6. ELASTIC SOLVER (Boundary Conditions: Fix left, Pull right)
    pull_force = 0.05 # Applied force in Newtons
    println("Solving elastic system (Pull Force: $(pull_force) N)...")
    
    # The solver handles Dirichlet boundaries directly on the sparse matrix
    U = solve_elastic_system(K_global, nodes, pull_force)
    
    # 7. EXTRACT DEFORMED GEOMETRY
    # The displacement vector U contains 3 DOFs per node: (u_x, u_y, theta)
    deformed_nodes = Vector{Tuple{Float64, Float64}}(undef, length(nodes))
    magnification_factor = 20.0 # Amplify deformations purely for visual clarity
    
    max_disp = 0.0
    for i in 1:length(nodes)
        u_x = U[3*i - 2]
        u_y = U[3*i - 1]
        
        # Calculate maximum true displacement for clinical logging
        max_disp = max(max_disp, sqrt(u_x^2 + u_y^2))
        
        # Apply amplified displacements to original coordinates
        deformed_nodes[i] = (nodes[i][1] + u_x * magnification_factor, 
                             nodes[i][2] + u_y * magnification_factor)
    end
    
    println("Maximum true displacement: ", round(max_disp * 1000, digits=4), " mm")

    # 8. VISUALIZATION (Overlay Original vs Deformed)
    println("\nRendering expansion plot...")
    p = plot(title="Linear Elastic Expansion (x$(magnification_factor) Mag)", 
             aspect_ratio=:equal, 
             legend=false, 
             grid=false,
             axis=false,
             background_color=:white,
             size=(900, 500))

    # Plot original undeformed mesh (Light gray, dashed)
    for (n1, n2) in elements
        p1, p2 = nodes[n1], nodes[n2]
        plot!(p, [p1[1], p2[1]], [p1[2], p2[2]], 
              color=:lightgray, linewidth=1.0, linestyle=:dash)
    end

    # Plot deformed mesh (Crimson red, solid)
    for (n1, n2) in elements
        p1, p2 = deformed_nodes[n1], deformed_nodes[n2]
        plot!(p, [p1[1], p2[1]], [p1[2], p2[2]], 
              color=:crimson, linewidth=2.0)
    end

    # Save the output image securely to the results folder
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    filepath = joinpath(results_dir, "linear_expansion.png")
    savefig(p, filepath)
    
    println("Plot saved successfully to: ", filepath)
    println("Done!")
    
    # Explicitly return nothing to prevent VS Code interception
    return nothing
end

# Execute the main function
run_linear_expansion()