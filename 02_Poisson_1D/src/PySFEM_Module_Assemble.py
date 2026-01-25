import numpy as np

from PyFEM_Module_GaussShape import *
from PyFEM_Module_SourceExact import *

'''
================================================================
            ASSEMBLE SYSTEM AND APPLY DIRICHLET BCs
================================================================
'''
def assem_SFEM(u, ndim, nel, ndf, lint, num_elements, num_nodes, coords, conn, boundary_nodes_id, eps):

    KMatrix = np.zeros((num_nodes * ndf, num_nodes * ndf), dtype=np.float64)
    FVector = np.zeros(num_nodes * ndf, dtype=np.float64)

    # Gauss points and weights for integration
    cg, wg = gauss(lint)

    for elem in range(num_elements):
        #print('Element:', elem)
        nodes_on_elem = conn[elem]
        xl = np.zeros((ndim, nel), dtype=np.float64)
        for i in range(nel):
            xl[0,i] = coords[nodes_on_elem][i]
        #print(xl)

        Ke = np.zeros((ndf * nel, ndf * nel), dtype=np.float64)
        Fe = np.zeros(ndf * nel, dtype=np.float64)

        for q in range(lint):
            xi= cg[q]
            if nel == 2:
                shp, detjl = shape_lin(xi, xl)
            elif nel == 3:
                shp, detjl = shape_quad(xi, xl)

            dvol = detjl * wg[q]

            # Compute gradient of u at this quadrature point
            xq = 0.0
            uh = 0.0
            gradu = np.zeros(ndim, dtype=np.float64)
            for i in range(nel):
                uh += shp[0, i] * u[nodes_on_elem[i]]
                gradu[0] += shp[1, i] * u[nodes_on_elem[i]]
                xq += shp[0,i] * xl[0,i]

            force_func = source(xq, eps)

            # Assemble stiffness matrix and residual vector correctly:
            for i in range(nel):
                for j in range(nel):
                    Ke[i * ndf, j * ndf] += shp[1, i] * shp[1, j] * dvol 

            for i in range(nel):
                # Residual form for nonlinear problem is typically: R_i = \int grad(phi_i) \cdot grad(u) dx - f_i
                # Here, approximating residual vector Fe (negative gradient of energy functional)
                Fe[i * ndf] += -gradu[0] * shp[1, i] * dvol + shp[0,i] * force_func * dvol 

        # Map local to global DOFs
        dof_map = []
        for n in nodes_on_elem:
            for d in range(ndf):
                dof_map.append(n * ndf + d)

        # Assemble into global KMatrix and FVector
        for a_local in range(ndf * nel):
            a_global = dof_map[a_local]
            FVector[a_global] += Fe[a_local]
            for b_local in range(ndf * nel):
                b_global = dof_map[b_local]
                KMatrix[a_global, b_global] += Ke[a_local, b_local]

    # Apply Dirichlet BC strongly in both matrix and residual
    for node_id in boundary_nodes_id:
        dof_idx = node_id * ndf
        KMatrix[dof_idx, :] = 0.0
        KMatrix[:, dof_idx] = 0.0
        KMatrix[dof_idx, dof_idx] = 1.0
        FVector[dof_idx] = exact_solution(coords[node_id], eps)  # Set exact BC value here

    return KMatrix, FVector

