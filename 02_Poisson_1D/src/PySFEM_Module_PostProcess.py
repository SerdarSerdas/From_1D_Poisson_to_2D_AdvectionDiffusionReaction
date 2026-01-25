import numpy as np
from PyFEM_Module_GaussShape import *
from PyFEM_Module_SourceExact import *
import csv

def compute_L2error(lint, ndim, nel, ndf, fem_method, InterpolOrder, mesh, eps):

    l2_error_u_squared = 0.0
    l2_error_q_squared = 0.0

    # Load CSV assuming numeric data, with a header line
    if InterpolOrder == 'P1':
        data = np.genfromtxt(f'Solution_{fem_method}_{InterpolOrder}_{mesh}.csv', delimiter=',', dtype=[('id_element', int), ('id_node1', int),  ('id_node2', int), ('x_node1', np.float64), ('x_node2', np.float64), ('u_node1', np.float64), ('u_node2', np.float64)], skip_header=1)
    
        num_elements = data.shape[0]

        output_data = []

        cg, wg = gauss(lint)

        for elem in range(data.shape[0]):
            #print(data['x_node1'][elem], data['x_node2'][elem])
            xl = np.zeros((ndim, nel), dtype=np.float64)
            ul = np.zeros((ndim, nel), dtype=np.float64)
            xl[0,0] = data['x_node1'][elem]
            xl[0,1] = data['x_node2'][elem]
            ul[0,0] = data['u_node1'][elem]
            ul[0,1] = data['u_node2'][elem]


            for q in range(lint):
                xi= cg[q]
                shp, detjl = shape_lin(xi, xl)
                dvol = detjl * wg[q]

                xq = 0.0
                uq = 0.0
                graduq = 0.0
                for i in range(nel):
                    xq     += shp[0, i] * xl[0,i]
                    uq     += shp[0, i] * ul[0,i]
                    graduq += shp[1, i] * ul[0,i]

                # output of flux
                if mesh == 1024 and q == 2:
                    print(xq, graduq)
                    output_data.append([xq, graduq])
                
                ue_q = exact_solution(xq, eps)
                qe_q = exact_deriv(xq, eps)
                diff_u = (ue_q - uq)**2 * dvol
                diff_q = (qe_q - graduq)**2 * dvol

                l2_error_u_squared += diff_u
                l2_error_q_squared += diff_q

        if mesh == 1024:
            with open(f'Gradient_data_{fem_method}_{InterpolOrder}_{mesh}.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['x_q', 'graduq'])  # Header
                writer.writerows(output_data)
    
    
    elif InterpolOrder == 'P2':
        data = np.genfromtxt(f'Solution_{fem_method}_{InterpolOrder}_{mesh}.csv', delimiter=',', dtype=[('id_element', int), ('id_node1', int),  ('id_node2', int),  ('id_node3', int), ('x_node1', np.float64), ('x_node2', np.float64), ('x_node3', np.float64), ('u_node1', np.float64), ('u_node2', np.float64), ('u_node3', np.float64)], skip_header=1)

        num_elements = data.shape[0]

        output_data = []

        cg, wg = gauss(lint)

        for elem in range(data.shape[0]):
            #print(data['x_node1'][elem], data['x_node2'][elem])
            xl = np.zeros((ndim, nel), dtype=np.float64)
            ul = np.zeros((ndim, nel), dtype=np.float64)
            xl[0,0] = data['x_node1'][elem]
            xl[0,1] = data['x_node2'][elem]
            xl[0,2] = data['x_node3'][elem]
            ul[0,0] = data['u_node1'][elem]
            ul[0,1] = data['u_node2'][elem]
            ul[0,2] = data['u_node3'][elem]

            for q in range(lint):
                xi= cg[q]
                shp, detjl = shape_quad(xi, xl)

                dvol = detjl * wg[q]

                xq = 0.0
                uq = 0.0
                graduq = 0.0
                for i in range(nel):
                    xq     += shp[0, i] * xl[0,i]
                    uq     += shp[0, i] * ul[0,i]
                    graduq += shp[1, i] * ul[0,i]

                # output of flux
                if mesh == 1024 and q == 2:
                    print(xq, graduq)
                    output_data.append([xq, graduq])
                
                ue_q = exact_solution(xq, eps)
                qe_q = exact_deriv(xq, eps)
                diff_u = (ue_q - uq)**2 * dvol
                diff_q = (qe_q - graduq)**2 * dvol

                l2_error_u_squared += diff_u
                l2_error_q_squared += diff_q

        if mesh == 1024:
            with open(f'Gradient_data_{fem_method}_{InterpolOrder}_{mesh}.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['x_q', 'graduq'])  # Header
                writer.writerows(output_data)

    print(f"L2 error u for {num_elements} elements:      {np.sqrt(l2_error_u_squared):.4e}")
    print(f"L2 error q for {num_elements} elements:      {np.sqrt(l2_error_q_squared):.4e}")

    return np.sqrt(l2_error_u_squared), np.sqrt(l2_error_q_squared)