import numpy as np
from matplotlib import pyplot as plt
from PyFEM_Module_SourceExact import *

def vizmesh(bc_coord_left, bc_coord_right, coords, nodes, List_boundary_id):
    y = np.linspace(bc_coord_left, bc_coord_right, 1)
    X, Y = np.meshgrid(coords, y)
    plt.plot(X, Y, 'k.', markersize=3)  # plot grid points as black dots
    plt.xlabel('x')
    plt.ylabel('y')
    # Mark the first and last x-coordinate with different markers
    plt.plot(nodes['x_node'][List_boundary_id[0]], y[0], 'ro', label='First point')   # red circle marker
    plt.plot(nodes['x_node'][List_boundary_id[1]], y[-1], 'ro', label='Last point')  # green square marker
    plt.title('Mesh grid of coords and y')
    plt.show()
    return

def plot_result(coords, u, num_elements, InterpolOrder, eps):
    xe = np.linspace(0,1,4001)
    plt.rcParams.update({'font.size': 14})
    plt.rcParams["figure.figsize"] = [14,6]
    plt.subplot(1,1,1)
    Result_x_u_unsorted = np.column_stack((coords, u))
    Result_x_u_sorted = Result_x_u_unsorted[np.argsort(Result_x_u_unsorted[:, 0])]
    plt.plot(xe, exact_solution(xe, eps), '-', color = 'red',   linewidth=2.5, label=f'Exact')
    plt.plot(Result_x_u_sorted[:,0], Result_x_u_sorted[:,1], '-o',   color = 'black', markersize = 3., linewidth=0.75, label=f'SFEM {InterpolOrder}')
    plt.xlabel('x')
    plt.ylabel(r'u(x)')
    plt.title(r'n$_\text{el}$=' f'{num_elements}')
    plt.legend()
    plt.grid(True)
    plt.show()
    return