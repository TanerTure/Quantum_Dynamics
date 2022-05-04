#Implementation of matrix numerov method, as described in Pillai2012

import numpy as np

def matrix_numerov(x_vals, potential_vals, order=1):
    size = x_vals.size - 2
    A_coefficients, B_coefficients = get_numerov_coefficients(order=order)

    operator = calculate_operator(x_vals, potential_vals, 
                                A_coefficients, B_coefficients, size)
    
    evals, evecs = diagonalize_matrix(operator)
    return evals, evecs




def get_numerov_coefficients(order=1):
    if(order == 1):
        A_coefficients = np.array([1, -2, 1], dtype=np.complex128)
        B_coefficients = 1/12 * np.array([1, 10, 1], dtype=np.complex128)
    
    if(order == 2):
        A_coefficients = -1/12 * np.array([1, -16, 30, -16, 1], dtype=np.complex128)
        B_coefficients =1/90 * np.array([1, -4, 96, -4, 1], dtype=np.complex128)
  
    if(order == 3):
        A_coefficients = 1/180 * np.array([2, -27, 270, -490,
                                        270, -27, 2], dtype=np.complex128)
        B_coefficients = 1/560 * np.array([1, -6, 15, 540,
                                        15, -6, 1], dtype=np.complex128)
    return (A_coefficients, B_coefficients)


def make_matrix(coefficients, size):
    matrix = np.zeros([size, size], dtype=np.complex128)
    diag = np.arange(size)
    num_diagonals = coefficients.size
    
    if(num_diagonals % 2 == 1):
        middle_index =int((coefficients.size-1)/2)
        np.fill_diagonal(matrix, coefficients[middle_index])
   
        for i in range(middle_index):
            offset = np.abs(i - middle_index) 
            matrix[diag[:-offset]+offset, diag[:-offset]] = coefficients[i]
            matrix[diag[:-offset], diag[:-offset]+offset] = coefficients[-i-1]
    return matrix

def calculate_operator(x_vals, potential_vals, A_coefficients, B_coefficients, size):
    A_matrix = make_matrix(A_coefficients, size)
    B_matrix = make_matrix(B_coefficients, size)
    
    stepsize = x_vals[1] - x_vals[0]
    B_matrix *= -2 * stepsize**2 #or -2 * m if mass is not 1 
    B_matrix = np.linalg.inv(B_matrix)

  
    diag = np.arange(size)
    V_matrix = np.zeros([size,size], dtype=np.complex128)
    np.fill_diagonal(V_matrix, potential_vals[1:-1])
   
    operator = (B_matrix @ A_matrix) + V_matrix

    return operator

def diagonalize_matrix(operator, sort=True):
    evals, evecs = np.linalg.eig(operator)
    if(sort == False):
        return evals, evecs
    else:
        energies = np.argsort(evals)
        evals = evals[energies]
        evecs[:,:] = evecs[:, energies]

        return evals, evecs

    
    
