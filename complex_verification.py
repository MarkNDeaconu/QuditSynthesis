import numpy as np
from qu5it import mat

zeta = np.exp(2j * np.pi / 5)
H_comp = np.array([[zeta**(i * j) for j in range(5)] for i in range(5)]) / np.sqrt(5)

T_comp = np.diag([1, zeta, zeta**3, zeta**2, zeta**4])

result = np.dot(H_comp, T_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, T_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, T_comp)
result = np.dot(result, T_comp)
result = np.dot(result, H_comp)
result = np.dot(result, T_comp)
result = np.dot(result, T_comp)

def are_matrices_equal(matrix1, matrix2, rtol=1e-05, atol=1e-08):
    return np.allclose(matrix1, matrix2, rtol=rtol, atol=atol)

print(are_matrices_equal(result, mat.comp()))
#print(H_comp*T_comp*T_comp*H*T_comp*H*T_comp*H*T_comp*T_comp*H*H*T_comp*H*T_comp*H*T_comp*T_comp*H*T_comp*H*T*T*T*H*T*T*T*H*T*T)
