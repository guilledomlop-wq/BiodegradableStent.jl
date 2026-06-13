"""
Example 01: Auxetic Stent Topology Generation

This script demonstrates how to generate and visualize the base 2D parametric mesh 
for the biodegradable stent. It highlights the use of the MeshGenerator module
and calculates standard clinical metrics. All visual results are saved to the 
local 'results' directory.
"""

# 1. ACTIVATE AND INSTANTIATE THE LOCAL ENVIRONMENT
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# 2. LOAD THE PACKAGE NATIVELY
# Thanks to Project.toml, we can load this exactly like a registered Julia package.
using BiodegradableStent
using Plots

# Force Plots (GR backend) into headless mode to prevent any window pop-ups in VS Code/REPL
ENV["GKSwstype"] = "100"

function run_topology_example()
    println("Starting Auxetic Topology Generation...")

    # 3. PARAMETRIC DEFINITION
    # Grid dimensions
    rows = 4
    cols = 8
    
    # Unit cell geometry
    h = 0.002         # Vertical strut height (m)
    l = 0.002         # Diagonal strut length (m)
    theta_deg = -15.0 # Re-entrant angle for auxetic behavior (negative for auxetic)
    
    # Cross-sectional dimensions
    width = 0.0001    # In-plane thickness (m)
    depth = 0.0001    # Out-of-plane thickness (m)

    # 4. PROPERTY CALCULATION
    A, I = get_cross_section_props(width, depth)
    println("Calculated Cross-Sectional Area: ", A, " m^2")
    println("Calculated Area Moment of Inertia: ", I, " m^4")

    # 5. MESH GENERATION
    println("Generating staggered mesh for $(rows)x$(cols) cells...")
    nodes, elements, w_array, A_array, I_array = generate_mesh_topology(
        rows, cols, h, l, theta_deg, width, A, I
    )

    println("Mesh generated successfully:")
    println("  -> Total Nodes: ", length(nodes))
    println("  -> Total Struts (Elements): ", length(elements))

    # 6. CLINICAL METRICS
    # Calculate total mass assuming the WE43 Magnesium alloy
    vol_m3, mass_mg = calculate_stent_mass(nodes, elements, A_array, WE43.rho)
    println("\nStent Clinical Metrics (Material: ", WE43.name, "):")
    println("  -> Volume: ", round(vol_m3 * 1e9, digits=3), " mm^3")
    println("  -> Mass:   ", round(mass_mg, digits=2), " mg")

    # 7. VISUALIZATION
    println("\nRendering plot...")
    
    # Initialize an empty plot with clinical aesthetic settings
    p = plot(title="Auxetic Stent Topology (2D Unrolled)", 
             aspect_ratio=:equal, 
             legend=false, 
             grid=false,
             axis=false,
             background_color=:white,
             size=(800, 400))

    # Plot all structural struts (Elements)
    for (n1, n2) in elements
        p1, p2 = nodes[n1], nodes[n2]
        plot!(p, [p1[1], p2[1]], [p1[2], p2[2]], 
              color=:steelblue, linewidth=1.5, alpha=0.9)
    end

    # Plot junctions (Nodes)
    x_coords = [n[1] for n in nodes]
    y_coords = [n[2] for n in nodes]
    scatter!(p, x_coords, y_coords, 
             color=:darkred, markersize=2.5, markerstrokewidth=0)

    # Save the output image securely to the results folder
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    filepath = joinpath(results_dir, "auxetic_mesh.png")
    savefig(p, filepath)
    
    println("Plot saved successfully to: ", filepath)
    println("Done!")
    
    # Explicitly return nothing so VS Code does not catch and display the plot object
    return nothing
end

# Execute the main function
run_topology_example()