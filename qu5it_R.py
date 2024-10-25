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




with open('cliffords5.pkl', 'rb') as f:
    cliffords = pickle.load(f)

full_set = [a * b  for a in [I,B,D] for b in cliffords]
        
edges = [I,B,D]


reduced_cyclotomics = []

for x in range(5):
    for y in range(5):
        for z in range(5):
            for w in  range(5):
                cyc = cyclotomic_element(z5, [x, y, z, w, 0],10)

                if cyc.sde == 10:
                    reduced_cyclotomics.append(cyc)


print(len(reduced_cyclotomics))

for i in range(3):
    mat= from_orbit([H,S,T])
    print(mat.pmap_state())
    for sta in mat.pmap_state():
        print(synth_search(sta, full_set))

    print(mat.sde_profile())


# for a in reduced_cyclotomics:
#     for b in reduced_cyclotomics:
#         for c in reduced_cyclotomics:
#             for d in reduced_cyclotomics:
#                 for e in reduced_cyclotomics:
#                     state = operator(5,1,[[a],[b],[c],[d],[e]])

#                     try:
#                         synth_search(state, full_set)
#                     except Exception:
#                         print('fail')
#             print('hi')


# print(len([operator(5,1,[[e0],[e0],[e0],[d],[e]]) for d in reduced_cyclotomics for e in reduced_cyclotomics]))

#eliminate clifford dupes

# eliminate state multiples

# for i in range(100):
#     orbit= set()

#     a = random.choice(reduced_cyclotomics)
#     b = random.choice(reduced_cyclotomics)
#     c = random.choice(reduced_cyclotomics)
#     d = random.choice(reduced_cyclotomics)
#     e = random.choice(reduced_cyclotomics)
#     sta = operator(5,1,[[a],[b],[c],[d],[e]])


#     if synth_search(sta, cliffords) != None: 
#         print('droppable')

#     else:
#         for cliff in cliffords:
#             new_state = (cliff* sta).pmap_state()[0]
#             if new_state.sde == 10:
#                 orbit.add(new_state)

#         print(len(orbit))

# for res in all_states:
#     try:
#         synth_search(res, full_set)
#     except Exception:
#         print('fail')

# cyc =  reduced_cyclotomics[10]

# print(cyc)