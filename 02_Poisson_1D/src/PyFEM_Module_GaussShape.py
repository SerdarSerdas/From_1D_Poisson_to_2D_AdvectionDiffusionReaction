import numpy as np


'''
========================================================
        Gauss Points/Weights & Ansatz Functions
========================================================
'''
def gauss(lint):
    if lint == 1:
        cg = np.array([0.0], dtype = float)
        wg = np.array([2.0], dtype = float)
    elif lint == 2:
        cg = np.array([-1./np.sqrt(3.0), 1./np.sqrt(3.0)], dtype = float)
        wg = np.array([1.0, 1.0], dtype = float)
    elif lint == 3:
        cg = np.array([-np.sqrt(3/5), 0.0, np.sqrt(3/5)], dtype = float)
        wg = np.array([5/9, 8/9, 5/9], dtype  = float) 
    elif lint == 4:
        cg = np.array([-0.8611363115940526, 
                       -0.3399810435848563, 
                       0.3399810435848563, 
                       0.8611363115940526], dtype = float)
        wg = np.array([0.3478548451374538, 
                       0.6521451548625461, 
                       0.6521451548625461, 
                       0.3478548451374538], dtype = float)
    return cg, wg

def shape_lin(xi, xl):
    shpl = np.zeros((2,2), dtype = float)
    shpl[0,0] = 0.5*(1.0 - xi)
    shpl[0,1] = 0.5*(1.0 + xi)
    sh0l = np.array([-0.5, 0.5], dtype = float)
    xjl = sh0l[0] * xl[0,0] + sh0l[1] * xl[0,1]
    detjl = xjl
    xjinv = 1.0/xjl
    shpl[1,0] = xjinv*sh0l[0]
    shpl[1,1] = xjinv*sh0l[1]
    return shpl, detjl

def shape_quad(xi, xl):
    shpq = np.zeros((2,3))
    shpq[0,0] = 0.5 * xi * (xi - 1.0)   # node 1
    shpq[0,1] = 0.5 * xi * (xi + 1.0)   # node 2
    shpq[0,2] = 1.0 - xi**2.0           # node 3
    sh0q = np.array([xi-0.5, xi+0.5, -2.0*xi])
    xjq = sh0q[0] * xl[0,0] + sh0q[1] * xl[0,1]  + sh0q[2] * xl[0,2]
    detjq = xjq
    xjinv = 1.0/xjq
    shpq[1,0] = xjinv*sh0q[0]
    shpq[1,1] = xjinv*sh0q[1]
    shpq[1,2] = xjinv*sh0q[2]
    return shpq, detjq
