import numpy as np
import pyamg
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import cg
from scipy.sparse import diags
from scipy.sparse.linalg import LinearOperator, cg


def solver_direct(KMatrix_sparse, FVector):
    '''
    Direct solver
    '''
    u = spsolve(KMatrix_sparse, FVector)
    residual = FVector - KMatrix_sparse @ u          # Compute residual
    residual_norm = np.linalg.norm(residual)  # Compute L2-norm Residual
    print(f"Residual norm after direct solve: {residual_norm:.6e}")
    return u


num_iterations = 0
residuals_cg = []


def solver_iterative(KMatrix_sparse, FVector, solver_type):
    '''
    Iterative solvers: CG without Pre-Conditioning
    '''
    num_iterations = 0
    residuals_cg = []

    def iteration_counter(xk):
        #global num_iterations
        nonlocal num_iterations, residuals_cg  # Reference LOCAL variables, NOT global
        num_iterations += 1
        r = FVector - KMatrix_sparse @ xk
        res_norm = np.linalg.norm(r)
        residuals_cg.append(res_norm)
        print(f"Iteration number {num_iterations}: Residual norm = {res_norm:8e}")

    x0 = np.zeros_like(FVector)

    # CG-method without pre-conditioning
    if solver_type == 'CG':     
        u, info = cg(KMatrix_sparse, FVector, x0=x0, rtol=1e-10, maxiter=100000, callback=iteration_counter)

    # CG-method with JACOBI pre-conditioning
    elif solver_type == 'CG-JACOBI': 
        # Construct Jacobi preconditioner with inverse of diagonal elements
        M_diag = 1.0 / KMatrix_sparse.diagonal()
        M = diags(M_diag)
        # Define a linear operator for the preconditioner
        M_inv = LinearOperator(KMatrix_sparse.shape, matvec=lambda x: M @ x)
        u, info = cg(KMatrix_sparse, FVector, M=M_inv, x0=x0, rtol=1e-10, maxiter=100000, callback=iteration_counter)

    elif solver_type == 'AMG-PCG':
        # Customize pre- and post-smoothers
        #presmoother = ('jacobi', {'omega': 2/3, 'iterations': 4})
        #postsmoother = ('gauss_seidel', {'sweep': 'symmetric', 'iterations': 4}) # 2707 / 593

        presmoother = ('gauss_seidel', {'sweep': 'symmetric', 'iterations': 4})
        postsmoother = ('gauss_seidel', {'sweep': 'symmetric', 'iterations': 4}) # 321 / 143

        #presmoother = ('jacobi', {'omega': 0.7, 'iterations': 6})
        #postsmoother = ('jacobi', {'omega': 0.7, 'iterations': 6}) # 613  / 284

        # Create an AMG preconditioner for the sparse matrix
        #ml = pyamg.ruge_stuben_solver(KMatrix_sparse, presmoother=presmoother, postsmoother=postsmoother)
        ml = pyamg.smoothed_aggregation_solver(KMatrix_sparse, presmoother=('gauss_seidel', {'sweep': 'symmetric', 'iterations': 4}), postsmoother=('gauss_seidel', {'sweep': 'symmetric', 'iterations': 4}))
        M = ml.aspreconditioner(cycle='W')
        # Solve with outer CG using AMG as preconditioner
        u, info = cg(KMatrix_sparse, FVector, M=M, tol=1e-10, maxiter=100000, callback=iteration_counter)




    if info == 0:
        print("CG solver converged successfully.")
    elif info > 0:
        print(f"CG solver did not converge after {info} iterations.")
    else:
        print("CG solver failed.")
    
    print(f"CG iterations count: {num_iterations}")

    return u, num_iterations