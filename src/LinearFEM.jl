module LinearFEM

export build_local_stiffness, assemble_global_matrix, solve_elastic_system

using LinearAlgebra
using SparseArrays 

# 1. LOCAL STIFFNESS CALCULATION
"""
Calculates the local stiffness matrix for a single 2D beam element.
- E: Young's Modulus
- A: Cross-sectional area
- I: Area moment of inertia
- L: Length of the strut
"""
function build_local_stiffness(E::T, A::T, I::T, L::T) where T<:Real
    c1 = (E * A) / L
    c2 = (T(12.0) * E * I) / (L^3)
    c3 = (T(6.0) * E * I) / (L^2)
    c4 = (T(4.0) * E * I) / L
    c5 = (T(2.0) * E * I) / L

    k_local = [
          c1    T(0.0)  T(0.0)  -c1    T(0.0)  T(0.0);
          T(0.0)  c2     c3    T(0.0) -c2     c3;
          T(0.0)  c3     c4    T(0.0) -c3     c5;
         -c1    T(0.0)  T(0.0)   c1    T(0.0)  T(0.0);
          T(0.0) -c2    -c3    T(0.0)  c2    -c3;
          T(0.0)  c3     c5    T(0.0) -c3     c4
    ]
    
    return k_local
end

# 2. GLOBAL MATRIX ASSEMBLY (HPC Optimized)
"""
Assembles the global stiffness matrix using SparseArrays for extreme memory efficiency.
Avoids dense matrix allocation to handle highly refined stent topologies.
"""
function assemble_global_matrix(nodes, elements, E_array::Vector{T}, A_array::Vector{T}, I_array::Vector{T}) where T<:Real
    total_dof = length(nodes) * 3 
    num_elements = length(elements)
    
    expected_entries = num_elements * 36
    I_idx = Vector{Int}(undef, expected_entries)
    J_idx = Vector{Int}(undef, expected_entries)
    V_val = Vector{T}(undef, expected_entries)
    
    pointer = 1
    
    for (i, (n1, n2)) in enumerate(elements)
        E_effective = E_array[i]
        A_loc = A_array[i]
        I_loc = I_array[i]

        p1, p2 = nodes[n1], nodes[n2]
        dx, dy = p2[1] - p1[1], p2[2] - p1[2]
        L = sqrt(dx^2 + dy^2)
        c, s = dx / L, dy / L 
        
        c1 = (E_effective * A_loc) / L
        c2 = (T(12.0) * E_effective * I_loc) / (L^3)
        c3 = (T(6.0) * E_effective * I_loc) / (L^2)
        c4 = (T(4.0) * E_effective * I_loc) / L
        c5 = (T(2.0) * E_effective * I_loc) / L
         
        cc, ss, cs = c*c, s*s, c*s

        k11 = c1*cc + c2*ss
        k12 = (c1 - c2)*cs
        k13 = -c3*s
        k14 = -k11
        k15 = -k12
        k16 = k13

        k22 = c1*ss + c2*cc
        k23 = c3*c
        k24 = -k12
        k25 = -k22
        k26 = k23

        k33 = c4
        k34 = c3*s
        k35 = -c3*c
        k36 = c5

        k44 = k11
        k45 = k12
        k46 = -k13

        k55 = k22
        k56 = -k23

        k66 = c4

        k_global = (
            k11, k12, k13, k14, k15, k16,
            k12, k22, k23, k24, k25, k26,
            k13, k23, k33, k34, k35, k36,
            k14, k24, k34, k44, k45, k46,
            k15, k25, k35, k45, k55, k56,
            k16, k26, k36, k46, k56, k66
        )
        
        dofs = (3*n1-2, 3*n1-1, 3*n1, 3*n2-2, 3*n2-1, 3*n2)
        
        idx = 1
        for row in 1:6
            for col in 1:6
                I_idx[pointer] = dofs[row]
                J_idx[pointer] = dofs[col]
                V_val[pointer] = k_global[idx]
                pointer += 1
                idx += 1
            end
        end
    end
    
    K = sparse(I_idx, J_idx, V_val, total_dof, total_dof)
    return K
end

# 3. ELASTIC SYSTEM SOLVER
"""
Applies boundary conditions and solves the system F = K * U for the elastic regime.
- K_global: The assembled sparse stiffness matrix.
- nodes: Array of (x,y) coordinates.
- pull_force: The total force applied to the right edge (Newtons).
"""
function solve_elastic_system(K_global, nodes, pull_force::T) where T<:Real
    total_dof = size(K_global, 1)
    F = zeros(T, total_dof)
    
    min_x, max_x = extrema(n -> n[1], nodes)
    
    # Boundary Band Tolerance (Fixes the "Hinge" effect on staggered grids)
    # Grabs the first and last 5% of the stent length to ensure the entire jagged edge is caught
    x_tol = (max_x - min_x) * T(0.05)
    
    left_nodes = findall(n -> n[1] <= min_x + x_tol, nodes)
    right_nodes = findall(n -> n[1] >= max_x - x_tol, nodes)
    
    force_per_node = pull_force / length(right_nodes)
    for n in right_nodes
        F[3*n - 2] = force_per_node 
    end
    
    is_fixed = falses(total_dof)
    for n in left_nodes
        is_fixed[3*n-2] = true
        is_fixed[3*n-1] = true
        is_fixed[3*n]   = true
    end
    
    for col in 1:total_dof
        col_fixed = is_fixed[col]
        for p in K_global.colptr[col]:(K_global.colptr[col+1]-1)
            row = K_global.rowval[p]
            if col_fixed || is_fixed[row]
                if row == col
                    K_global.nzval[p] = T(1.0)
                else
                    K_global.nzval[p] = T(0.0)
                end
            end
        end
    end
    
    for i in 1:total_dof
        if is_fixed[i]
            F[i] = T(0.0)
        end
    end
    
    return K_global \ F
end

end