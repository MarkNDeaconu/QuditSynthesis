from datastructures import *
import random
import pickle
import itertools
from multiprocessing import Pool


z5 =  cyclotomic_ring(5,math.sqrt(5))

n = cyclotomic_element(z5, [0,0,0,0,0]) 
e0 = cyclotomic_element(z5, [1,0,0,0,0])
e1 = cyclotomic_element(z5, [0,1,0,0,0])
e2 = cyclotomic_element(z5, [0,0,1,0,0])
e3 = cyclotomic_element(z5, [0,0,0,1,0])
e4 = cyclotomic_element(z5, [0,0,0,0,1])


H = (1/math.sqrt(5))*operator(5,5, [[e0,e0,e0,e0,e0], [e0,e1,e2,e3,e4], [e0,e2,e4,e1,e3], [e0,e3,e1,e4,e2], [e0, e4, e3,e2,e1]])
H.string = 'H'
R = operator(5,5,[
    [e0, n, n, n, n],  
    [n, e0, n, n, n], 
    [n, n, e0, n, n], 
    [n, n, n, e0, n], 
    [n, n, n, n, (-1)*e0],  
])


R.string = 'R'

T = operator(5,5, [
    [e0,n,n,n,n],
    [n,e1,n,n,n],
    [n,n,e3,n,n],
    [n,n,n,e2,n],
    [n,n,n,n,e4]
])

S= operator(5,5, [
    [e0,n,n,n,n],
    [n,e1,n,n,n],
    [n,n,e3,n,n],
    [n,n,n,e1,n],
    [n,n,n,n,e0]
])

S.string = 'S'

T.string = 'T'



def D_gate(a,b,c,d,e):
    power_map = [e0, e1, e2, e3, e4]
    gate = operator(5,5, [[power_map[a],n,n,n,n], [n,power_map[b],n,n,n], [n,n,power_map[c],n,n], [n,n,n,power_map[d],n], [n, n, n, n, power_map[e]]])
    gate.string = 'D('+str(a)+str(b)+str(c)+str(d)+str(e)+')'
    return(gate)

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R
I = R*R
I.string = ''




# with open('cliffords5.pkl', 'rb') as f:
#     cliffords = pickle.load(f)

# full_set = [a * b  for a in [I,B,D] for b in cliffords]
        
# edges = [I,B,D]

