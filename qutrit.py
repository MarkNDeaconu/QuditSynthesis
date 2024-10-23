from datastructures import *
import random
import pickle
import itertools


z3 =  cyclotomic_ring(3,complex(0, math.sqrt(3)))
n = cyclotomic_element(z3, [0,0,0]) 
e0 = cyclotomic_element(z3, [1,0,0])
e1 = cyclotomic_element(z3, [0,1,0])
e2 = cyclotomic_element(z3, [0,0,1])

e0m = cyclotomic_element(z3, [1,0,0],1)
e1m = cyclotomic_element(z3, [0,1,0],1)
e2m = cyclotomic_element(z3, [0,0,1],1)

H = operator(3,3, [
    [e0m,e0m,e0m], 
    [e0m,e1m,e2m], 
    [e0m,e2m,e1m]
])

R = operator(3,3,[
    [e0,n,n],
    [n,e0,n],
    [n,n,(-1)*e0]
])

S = operator(3,3,[
    [e0,n,n],
    [n,e1,n],
    [n,n,e0]
])


H.string = 'H'
R.string = 'R'
S.string = 'S'

I=R*R
I.string = ''

def D_gate(a,b,c):
    power_map = [e0, e1, e2]
    return(operator(3,3, [[power_map[a],n,n], [n,power_map[b],n], [n,n,power_map[c]]]  ))

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R

with open('cliffords3.pkl', 'rb') as f:
    cliffords = pickle.load(f)

full_set = [a * b  for a in [I,B] for b in cliffords]

edges = [I,B]

mat=from_orbit([H,S,R])

print(synth_search(mat, full_set))

print(neighbors_mat(mat, [I,B,D], full_set))


print(mat)

print(mat*B)
print(mat*D)





