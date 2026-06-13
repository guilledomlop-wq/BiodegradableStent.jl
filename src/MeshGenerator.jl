module MeshGenerator

export generate_mesh_topology, get_cross_section_props, calculate_stent_mass

# 1. MESH TOPOLOGY (Nodes and Elements)
"""
Generates the nodes and connectivity (elements) for the staggered auxetic mesh.
Merges duplicate nodes automatically to ensure a continuous solid structure.
- Returns:
  - nodes: Array of (x, y) coordinate tuples.
  - elements: Array of (node_i, node_j) index tuples representing the struts.
"""
function generate_mesh_topology(rows::Int, cols::Int, h::T, l::T, theta_deg::T, w::T=T(0.0), area::T=T(0.0), inertia::T=T(0.0)) where T<:Real
    # 1.1 MESH KINEMATICS: Define the "Stitch" spacing
    theta = deg2rad(theta_deg)
    dx = T(2.0) * l * cos(theta)
    dy = h + l * sin(theta)

    # Internal dimensions for node placement
    x_ext = l * cos(theta)
    y_int_sup = h/T(2.0) + l * sin(theta)
    y_int_inf = -h/T(2.0) - l * sin(theta)

    # 1.2 NODE MERGING STORY
    # To avoid "floating" struts, shared nodes between cells must have the same ID.
    # We use a Dictionary to map physical coordinates to a unique integer ID.
    nodes_dict = Dict{Tuple{T, T}, Int}()
    nodes = Vector{Tuple{T, T}}()
    sizehint!(nodes, rows * cols * 6) # OPTIMIZATION: Pre-allocate memory
    
    elements = Vector{Tuple{Int, Int}}()
    sizehint!(elements, rows * cols * 6)

    # Internal helper to assign and retrieve unique node IDs
    function get_node_id(x::T, y::T)::Int
        # Rounding is CRITICAL: prevents 1.000000000001 being different from 1.0.
        # We use 8 digits to handle micrometer precision in a meter-scale world.
        rx, ry = round(x, digits=8), round(y, digits=8)
        id = get(nodes_dict, (rx, ry), 0)
        if id == 0
            push!(nodes, (rx, ry))
            id = length(nodes)
            nodes_dict[(rx, ry)] = id
        end
        return id
    end

    # 1.3 GRID GENERATION LOOP
    for j in 1:rows
        for i in 1:cols
            # Offset logic for staggered rows (creates the "honeycomb" shift)
            offset_x = (j % 2 == 0) ? (dx / T(2.0)) : T(0.0)
            cx = (i - T(1.0)) * dx + offset_x
            cy = (j - T(1.0)) * dy

            # Get unique IDs for these 6 points of the current cell
            n1 = get_node_id(cx - x_ext, cy + h/T(2.0))
            n2 = get_node_id(cx - x_ext, cy - h/T(2.0))
            n3 = get_node_id(cx + x_ext, cy + h/T(2.0))
            n4 = get_node_id(cx + x_ext, cy - h/T(2.0))
            n5 = get_node_id(cx, cy + y_int_sup)
            n6 = get_node_id(cx, cy + y_int_inf)

            # 1.4 STRUT CONNECTIVITY
            # Using minmax ensures (1,2) and (2,1) are treated as the same strut.
            push!(elements, minmax(n1, n2)) # Left pillar
            push!(elements, minmax(n3, n4)) # Right pillar
            push!(elements, minmax(n1, n5)) # Top-left diagonal
            push!(elements, minmax(n3, n5)) # Top-right diagonal
            push!(elements, minmax(n2, n6)) # Bottom-left diagonal
            push!(elements, minmax(n4, n6)) # Bottom-right diagonal
        end
    end

    # Clean up directional duplicates and shared pillars between cells
    unique!(elements)
    
    num_elements = length(elements)
    w_array = fill(w, num_elements)
    area_array = fill(area, num_elements)
    inertia_array = fill(inertia, num_elements)

    return nodes, elements, w_array, area_array, inertia_array
end

# 2. CROSS-SECTIONAL PROPERTIES
"""
Calculates the Area (A) and Area Moment of Inertia (I) for the stent struts.
Assumes a rectangular cross-section.
- t: The in-plane strut thickness (width).
- depth: The out-of-plane thickness (stent wall).
"""
function get_cross_section_props(t::T, depth::T) where T<:Real
    A = t * depth # Resistance to axial pulling
    I = (depth * t^T(3.0)) / T(12.0) # Resistance to bending (core of stent collapse)
    return A, I
end

# 3. CLINICAL METRICS (Mass and Volume)
"""
Calculates the total volume and mass of the stent based on its topology.
- area_array: Vector of cross-sectional areas for each element (m^2).
- rho: Material density (kg/m^3).
Returns: (volume_m3, mass_mg).
"""
function calculate_stent_mass(nodes, elements, area_array::Vector{T}, rho::T) where T<:Real
    total_volume = T(0.0)
    # Logical Story: Iterate over all unique struts to sum up material volume.
    for (i, (n1, n2)) in enumerate(elements)
        p1, p2 = nodes[n1], nodes[n2]
        # Calculate Euclidean length L using hypot for numerical stability
        L = hypot(p2[1] - p1[1], p2[2] - p1[2])
        total_volume += L * area_array[i]
    end
    # Mass in milligrams (kg -> mg requires 1e6 factor)
    mass_mg = (total_volume * rho) * T(1e6)
    return total_volume, mass_mg
end

end