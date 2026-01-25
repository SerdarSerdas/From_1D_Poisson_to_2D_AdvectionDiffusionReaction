"""
SFEM:   STANDARD FEM (GALERKIN) 
LSFEM:  LEAST-SQUARES FEM

POISSON EQUATION --> 1D



01.11.25    Dr.-Ing. Serdar Serdas

The following naming convention of variables are used.


=====================================================================
Name                Description
=====================================================================

Input Parameters:
=================

fem_method          type = string:  chooses the method: 'SFEM' or 'LSFEM'
InterpolOrder       type = string:  'LINEAR'    --> nel = 2
                                    'QUADRATIC' --> nel = 3
solver_method       type = string:  'DIRECT' or 'ITERATIVE'

solver_type         'CG'           Pre-Cond: No
                    'CG-JACOBI'    Pre-Cond: Jacobi
                    'Iter-ILU'          Pre-Cond: Incomplete LU-Factorization
                    'Iter-AMG'          Pre-Cond: No
                    'Iter-AMG-Pre'      Pre-Cond: 

PostProc            True  = Generating results and save as PNG-File 
                    False = nothing will be done
bc_coord_left       coordinates of left boundary
bc_coord_right      coordinates of rightt boundary

----------------------------------------------------------------

Fixed Parameters:
=================

Linear Elements

SFEM:   ndm = 1     number of dimension
        ndf = 1     number of degree of freedom
        nel = 2     number of nodes per element

LSFEM:  ndm = 1     number of dimension
        ndf = 2     number of degree of freedom
        nel = 2     number of nodes per element

================================================================
"""

import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
from pathlib import Path
from scipy.sparse.linalg import eigsh

from PyFEM_Module_GaussShape import *
from PyFEM_Module_SourceExact import *
from PySFEM_Module_Assemble import *
from PyFEM_Module_Solvers import *
from PyFEM_Module_Vizualization import *
from PySFEM_Module_PostProcess import *

'''
================================================================
                    INPUT PARAMETERS
================================================================
'''
#Level        = [0, 1, 2, 3,  4,  5,  6,   7,   8,   9,   10]
#mesh_levels  = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

#mesh_levels   = [16, 32, 64, 128, 256]
mesh_levels     = [32, 64, 128, 256, 512, 1024]
eps             = 5.e-6          # Viscosity 1.e-3 or 1.e-5
InterpolOrder   = 'P1'          # 'P1' or 'P2'
grid            = 'RMR'
solver_method   = 'DIRECT'   # 'DIRECT' or 'ITERATIVE'
solver_type     = 'AMG-PCG'     # 'CG': CG-method, Pre-Conditioner:         NO
                                # 'CG-JACOBI': CG-method, Pre-Conditioner:  JACOBI
                                # 'AMG-PCG'
lint            = 3             # Number of Gausspoints
bc_coord_left   = 0.0           # Coordinate of left boundary
bc_coord_right  = 1.0           # Coordinate of right boundary
VizMesh         = False         # Visualize Mesh and Boundaries
PostProcess     = True         # For instance, computing L2-Error

MeshLevel_history    = []  # Mesh Level History
NumElement_history   = []  # Number of elements
NumNodes_history     = []  # Number of nodes
NumDof_history       = []  # Number of Dof
L2_error_u_history   = []  # L2-error for u
L2_error_q_history   = []  # L2-error for q 
Solver_iter_history  = []  # Iteration numbers of Solver
Condition_nr_history = []  # Condition number

for mesh in mesh_levels:
    '''
    ================================================================
                DEFINE THE PATH OF INPUT-FILE LOCATION
    ================================================================
    '''
    path_nodes = Path(f'C:/Users/serda/Desktop/GitHub/FEM/AbaqusMesh_1D_01/{grid}/PyNodes_{InterpolOrder}_{mesh}')
    path_elements = Path(f'C:/Users/serda/Desktop/GitHub/FEM/AbaqusMesh_1D_01/{grid}/PyElements_{InterpolOrder}_{mesh}')

    '''
    ================================================================
                FIXED PARAMETERS
    ================================================================
    '''
    fem_method      = 'SFEM'  # 'SFEM'
    ndim            = 1       # Number of dimension
    ndf             = 1

    if InterpolOrder == 'P1':
        nel = 2
        fmt = ('%d', '%d', '%d', '%.12e', '%.12e', '%.12e', '%.12e')
        output_header = 'ElementID,Node1ID,Node2ID,X_node1,X_node2,U_node1,U_node2'
    elif InterpolOrder == 'P2':
        nel = 3
        fmt = ('%d', '%d', '%d', '%d', '%.6e', '%.12e', '%.12e', '%.12e', '%.12e', '%.12e')
        output_header = 'ElementID,Node1ID,Node2ID,Node3ID,X_node1,X_node2,X_node3,U_node1,U_node2,U_node3'


    '''
    ================================================================
                REARRANGEMENT OF NODES & ELEMENTS
    ================================================================
    '''
    nodes    = np.loadtxt(path_nodes, dtype=[('id_node', int), ('x_node', np.float64)])
    elements = np.loadtxt(path_elements, dtype=int)

    coords = nodes['x_node']
    coords = coords.astype(np.float64)
    conn = elements[:, 1:]

    num_elements = elements.shape[0]  # total number of elements
    num_nodes    = nodes.shape[0]     # total number of nodes

    u = np.zeros((num_nodes*ndf), dtype=np.float64)

    '''
    ================================================================
                DEFINE BOUNDARY CONDITIONS
    ================================================================
    '''
    # Find indices of boundary nodes where coordinates are 0 and 1
    boundary_nodes_id = np.where((np.isclose(nodes['x_node'], bc_coord_left)) | (np.isclose(nodes['x_node'], bc_coord_right)))[0]
    boundary_coords = coords[boundary_nodes_id]
    boundary_coords.astype(np.float64)

    # Create a structured array for boundary data with integer IDs and float coords
    boundary_data = np.zeros(boundary_nodes_id.shape[0], dtype=[('id_node', int), ('x_node', np.float64)])
    boundary_data['id_node'] = nodes['id_node'][boundary_nodes_id]
    boundary_data['x_node'] = boundary_coords[:].astype(np.float64)

    # Assign exact solution values at boundary nodes as Dirichlet BC
    u[boundary_nodes_id] = exact_solution(nodes['x_node'][boundary_nodes_id], eps)
    u.astype(np.float64)

    #print(boundary_data['id_node'])
    #for idx in boundary_nodes_id:
    #    print(f"Node {idx}, x = {nodes['x_node'][idx]}, u = {u[idx]}")

    '''
    ================================================================
                VIZUALIZATION OF MESH & BOUNDARIES
    ================================================================
    '''
    List_boundary_id = []
    for idx in boundary_nodes_id:
        List_boundary_id.append(idx)

    if VizMesh:
        vizmesh(bc_coord_left, bc_coord_right, coords, nodes, List_boundary_id)


    '''
    ================================================================
                ASSEMBLE SYSTEM --> LHS & RHS
    ================================================================
    '''
    KMatrix, FVector = assem_SFEM(u, ndim, nel, ndf, lint, num_elements, num_nodes, coords, conn, boundary_nodes_id, eps)


    '''
    ================================================================
                SOLVER
    ================================================================
    '''

    KMatrix_sparse = csr_matrix(KMatrix)

    if solver_method == 'DIRECT':
        u = solver_direct(KMatrix_sparse, FVector)
        plot_result(coords, u, num_elements, InterpolOrder, eps)
    elif solver_method == 'ITERATIVE':  # Call CG without preconditioning (M=None)
        u, num_iterations = solver_iterative(KMatrix_sparse, FVector, solver_type)
        plot_result(coords, u, num_elements, InterpolOrder, eps)

    '''
    ================================================================
                OUTPUT SOLUTION-VECTOR U
    ================================================================
    '''
    output_data = []

    for elem_id in range(num_elements):
        node_ids = conn[elem_id]
        if InterpolOrder == 'P1':
            x1, x2 = coords[node_ids[0]], coords[node_ids[1]]
            u1, u2 = u[node_ids[0]], u[node_ids[1]]
            row = [elem_id, node_ids[0], node_ids[1], x1, x2, u1, u2]
        elif InterpolOrder == 'P2':
            x1, x2, x3 = coords[node_ids[0]], coords[node_ids[1]], coords[node_ids[2]]
            u1, u2, u3 = u[node_ids[0]], u[node_ids[1]], u[node_ids[2]]
            row = [elem_id, node_ids[0], node_ids[1], node_ids[2], x1, x2, x3, u1, u2, u3]
        output_data.append(row)

    # Convert to numpy array with dtype=object to hold mixed types safely
    output_array = np.array(output_data, dtype=object)
    np.savetxt(f'Solution_{fem_method}_{InterpolOrder}_{mesh}.csv', output_array, delimiter=',',header=output_header, comments='', fmt=fmt)

    print()
    print('='*50)
    print(f'SFEM - {InterpolOrder}')
    print('='*50)
    print(f'Mesh Level:                    {mesh}')
    print(f'Number of elements:            {num_elements}')
    print(f'Number of Nodes:               {num_nodes}')
    print(f'Number of Degree of Freedoms:  {num_nodes*ndf}')
    print()

    if PostProcess:
        error_u, error_q = compute_L2error(lint, ndim, nel, ndf, fem_method, InterpolOrder, mesh, eps)
        #print(f"L2 error u for {num_elements} elements:      {L2err_u:.4e}")
        #print(f"L2 error q for {num_elements} elements:      {L2err_q:.4e}")
    kappa = np.linalg.cond(KMatrix)
    print(f"Condition number (2-norm): {kappa:.1e}")

    # SAVE initial solve residuals
    MeshLevel_history.append(mesh)
    NumElement_history.append(num_elements)
    NumNodes_history.append(num_nodes)
    NumDof_history.append(num_nodes*ndf)
    L2_error_u_history.append(error_u)
    L2_error_q_history.append(error_q)
    Solver_iter_history.append(num_iterations)
    Condition_nr_history.append(kappa)

print("\n=== CONVERGENCE SUMMARY ===")
print("Mesh Level, #Elements, #Nodes, #DoF, L2-Error u, L2-Error q, iter_Solver, kappa")
for i in range(len(MeshLevel_history)):
    print(f"{MeshLevel_history[i]}, {NumElement_history[i]}, {NumNodes_history[i]}, {NumDof_history[i]}, {L2_error_u_history[i]:.4e}, {L2_error_q_history[i]:.4e}, {Solver_iter_history[i]}, {Condition_nr_history[i]:.1e}")

# NEW: Compute and print convergence rates
print("\n=== CONVERGENCE RATES ===")
print("Interval, u_ratio, u_order, q_ratio, q_order")
for i in range(1, len(L2_error_u_history)):
    u_ratio = L2_error_u_history[i-1] / L2_error_u_history[i]
    q_ratio = L2_error_q_history[i-1] / L2_error_q_history[i]
    u_order = np.log2(u_ratio)
    q_order = np.log2(q_ratio)
    print(f"{MeshLevel_history[i-1]}→{MeshLevel_history[i]}, {u_ratio:.2f}, {u_order:.2f}, {q_ratio:.2f}, {q_order:.2f}")

# Write convergence history to CSV
history_data = np.column_stack((
    np.array(MeshLevel_history),
    np.array(NumElement_history),
    np.array(NumNodes_history),
    np.array(NumDof_history),
    np.array(L2_error_u_history),
    np.array(L2_error_q_history),
    np.array(Solver_iter_history),
    np.array(Condition_nr_history)
))


np.savetxt(f'Conv_History_BL_SFEM_{grid}_{InterpolOrder}_{solver_type}.csv', 
           history_data, 
           delimiter=',',
           header='MeshLevel,NumElem,NumNodes,DoF,Error_u,Error_q,iter_Solver,kappa',
           comments='',
           fmt='%d,%d,%d,%d,%.4e,%.4e,%d,%.1e')

print("\nConvergence history saved to 'Convergence_History.csv'")