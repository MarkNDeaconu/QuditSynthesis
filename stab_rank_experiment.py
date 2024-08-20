from datastructures import *

from datastructures import *
import random
z8 =  cyclotomic_ring(8,math.sqrt(2))
n = cyclotomic_element(z8, [0,0,0,0,0,0,0,0]) 
e0 = cyclotomic_element(z8, [1,0,0,0,0,0,0,0])
e1 = cyclotomic_element(z8, [0,1,0,0,0,0,0,0])
e2 = cyclotomic_element(z8, [0,0,1,0,0,0,0,0])
e3 = cyclotomic_element(z8, [0,0,0,1,0,0,0,0])

T_eigen = state(2, [e0,(-1)*e0 + e1 + (-1)*e3])

print(T_eigen)

print(T_eigen.tensor_power(6))