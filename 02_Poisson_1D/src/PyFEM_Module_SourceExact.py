import numpy as np

'''
================================================================
Define the exact solution, the derivative and the source term
================================================================
'''

'''
EXACT SOLUTION
'''
def exact_solution(x,eps):
    return 1 - (np.sinh(x / np.sqrt(eps)) + np.sinh((1 - x) / np.sqrt(eps))) / np.sinh(1 / np.sqrt(eps))
    
'''
DERIVATIVE OF EXACT SOLUTION
'''
def exact_deriv(x,eps):
    u_prime = -(np.cosh(x / np.sqrt(eps)) / np.sqrt(eps) - np.cosh((1 - x) / np.sqrt(eps)) / np.sqrt(eps)) / np.sinh(1 / np.sqrt(eps))
    return u_prime

'''
SOURCE TERM
'''
def source(x,eps):
    u_primeprime = -(np.sinh(x / np.sqrt(eps)) / np.sqrt(eps) + np.sinh((1 - x) / np.sqrt(eps)) / np.sqrt(eps)) / (np.sinh(1 / np.sqrt(eps)) * np.sqrt(eps))
    
    return -u_primeprime
    
     
