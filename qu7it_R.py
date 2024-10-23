from datastructures import *
import random
import pickle


z7 =  cyclotomic_ring(7,complex(0, math.sqrt(7)))

n = cyclotomic_element(z7, [0,0,0,0,0,0,0]) 
e0 = cyclotomic_element(z7, [1,0,0,0,0,0,0])
e1 = cyclotomic_element(z7, [0,1,0,0,0,0,0])
e2 = cyclotomic_element(z7, [0,0,1,0,0,0,0])
e3 = cyclotomic_element(z7, [0,0,0,1,0,0,0])
e4 = cyclotomic_element(z7, [0,0,0,0,1,0,0])
e5 = cyclotomic_element(z7, [0,0,0,0,0,1,0])
e6 = cyclotomic_element(z7, [0,0,0,0,0,0,1])

e0m = cyclotomic_element(z7, [1,0,0,0,0,0,0],1)
e1m = cyclotomic_element(z7, [0,1,0,0,0,0,0],1)
e2m = cyclotomic_element(z7, [0,0,1,0,0,0,0],1)
e3m = cyclotomic_element(z7, [0,0,0,1,0,0,0],1)
e4m = cyclotomic_element(z7, [0,0,0,0,1,0,0],1)
e5m = cyclotomic_element(z7, [0,0,0,0,0,1,0],1)
e6m = cyclotomic_element(z7, [0,0,0,0,0,0,1],1)

H = operator(7,7, [
    [e0m, e0m, e0m, e0m, e0m, e0m, e0m],      
    [e0m, e1m, e2m, e3m, e4m, e5m, e6m],     
    [e0m, e2m, e4m, e6m, e1m, e3m, e5m],     
    [e0m, e3m, e6m, e2m, e5m, e1m, e4m],      
    [e0m, e4m, e1m, e5m, e2m, e6m, e3m],     
    [e0m, e5m, e3m, e1m, e6m, e4m, e2m],     
    [e0m, e6m, e5m, e4m, e3m, e2m, e1m]       
])

R = operator(7,7,[
    [e0, n, n, n, n, n, n],  
    [n, e0, n, n, n, n, n], 
    [n, n, e0, n, n, n, n], 
    [n, n, n, e0, n, n, n], 
    [n, n, n, n, e0, n, n],  
    [n, n, n, n, n, e0, n],  
    [n, n, n, n, n, n, (-1)*e0]  
])

I= R*R

S = operator(7,7,[
    [e0, n, n, n, n, n, n],  
    [n, e1, n, n, n, n, n], 
    [n, n, e3, n, n, n, n], 
    [n, n, n, e6, n, n, n], 
    [n, n, n, n, e3, n, n],  
    [n, n, n, n, n, e1, n],  
    [n, n, n, n, n, n, e0]  
])

H.string = 'H'
R.string = 'R'
S.string = 'S'

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R



with open('cliffords7fullset.pkl', 'rb') as f:
    full_set = pickle.load(f)

