"""
Example 03: Material Dispatch & Degradation Kinetics

This script demonstrates Julia's Multiple Dispatch architecture. 
It applies the Shrinking Core degradation model over a 6-month (180 days) 
clinical window, automatically routing the calculation to specific kinetic 
algorithms based on the abstract material type provided (Magnesium vs Zinc).
"""

# 1. ACTIVATE AND INSTANTIATE THE LOCAL ENVIRONMENT
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# 2. LOAD THE PACKAGE NATIVELY
using BiodegradableStent
using Plots

# Force Plots (GR backend) into headless mode to prevent window pop-ups
ENV["GKSwstype"] = "100"

function run_degradation_example()
    println("Starting Material Degradation & Dispatch Simulation...")

    # 3. CLINICAL TIMELINE & INITIAL STATE
    initial_width = 0.00015 # 150 microns starting thickness
    time_days = collect(0.0:5.0:180.0) # Evaluate every 5 days for 6 months
    delta_t = 5.0 # Time step in days

    # Standard corrosion rates for uniform degradation (mm/year)
    rate_WE43 = 0.020 # Magnesium degrades relatively fast
    rate_Zinc = 0.005 # Zinc has a much slower, controlled baseline rate

    # Arrays to store the history for plotting
    w_history_mg = Float64[initial_width]
    w_history_zn = Float64[initial_width]

    # Current geometry states (represented as arrays for the element loop)
    w_mg = [initial_width]
    w_zn = [initial_width]

    # 4. TIME INTEGRATION LOOP
    println("Simulating 180 days of chemical degradation...")
    
    for t in time_days[2:end]
        # =====================================================================
        # MULTIPLE DISPATCH IN ACTION:
        # The exact same function name is called, but Julia's compiler routes 
        # it to highly optimized, material-specific native code under the hood.
        # =====================================================================
        
        w_mg = calculate_corrosion_step(WE43, w_mg, delta_t, rate_WE43)
        w_zn = calculate_corrosion_step(PureZinc, w_zn, delta_t, rate_Zinc)
        
        # Log current thickness
        push!(w_history_mg, w_mg[1])
        push!(w_history_zn, w_zn[1])
    end

    println("\nFinal Clinical Status (Day 180):")
    println("  -> Thickness (", WE43.name, "): ", round(w_history_mg[end] * 1e6, digits=2), " µm")
    println("  -> Thickness (", PureZinc.name, "): ", round(w_history_zn[end] * 1e6, digits=2), " µm")

    # 5. VISUALIZATION
    println("\nRendering degradation kinetics plot...")
    
    p = plot(time_days, w_history_mg .* 1e6, 
             title="Material Kinetics (Multiple Dispatch)",
             xlabel="Time in-vivo (Days)", 
             ylabel="Strut Thickness (µm)",
             label=WE43.name, 
             color=:steelblue, lw=4, grid=true,
             size=(800, 500))

    plot!(p, time_days, w_history_zn .* 1e6, 
          label=PureZinc.name, 
          color=:darkgray, lw=4)

    # Highlight a hypothetical critical structural threshold
    hline!(p, [50.0], line=(2, :dash, :crimson), label="Critical Structural Threshold")

    # Save the output image securely to the results folder
    results_dir = joinpath(@__DIR__, "results")
    mkpath(results_dir)
    filepath = joinpath(results_dir, "material_degradation.png")
    savefig(p, filepath)
    
    println("Plot saved successfully to: ", filepath)
    println("Done!")
    
    return nothing
end

# Execute the main function
run_degradation_example()