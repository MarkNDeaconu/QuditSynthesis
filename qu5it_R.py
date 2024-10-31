from datastructures import *
import random
import pickle
import itertools
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor, as_completed

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

print((e2+e3 + (-1)*e4).norm())
# with open('cliffords5.pkl', 'rb') as f:
#     cliffords = pickle.load(f)

# full_set = [a * b  for a in [I,B,D] for b in cliffords]

# edges = [I,B,D]


# reduced_cyclotomics = []
# non_reduced = []



# diags = z5.torus(cliffords,n)

# full_set2 = [a * b * c for a in [I,B,D] for b in diags for c in [I,H, H*H, H*H*H]]


# for i in range(100):
#     mat = z5.from_orbit([H,S,R])

#     print(mat.synth_search(full_set2))

# for cliff in cliffords:
#     proper_rows = 0 
#     for row in cliff.matrix:
#         zeros = 0
#         for elem in row:
#             if elem == n:
#                 zeros+=1
#         if zeros == 4:
#             proper_rows+=1
#     if proper_rows ==5:
#         diags.append(cliff)


# print(len(diags))

        



# for w in  range(5):
#     cyc = cyclotomic_element(z5, [1, 1, 1, w, 0],10)
#     if cyc.sde == 10:
#         reduced_cyclotomics.append(cyc)

# possible_states = [operator(5,1,[[a],[b],[c],[d],[e]]) for a in reduced_cyclotomics for b in reduced_cyclotomics for c in reduced_cyclotomics for d in reduced_cyclotomics for e in reduced_cyclotomics]

# print(len(possible_states))

# def task(state):
#     if synth_search(state, full_set) == None:
#         return('fail')



# def run_parallel_task(states, max_workers=12):
#     with ProcessPoolExecutor(max_workers=max_workers) as executor:
#         futures = [executor.submit(task, state) for state in states]
        
#         for future in as_completed(futures):
#             try:
#                 future.result()  # Process result or handle exceptions
#             except Exception as e:
#                 print(f"Task generated an exception: {e}")

# print('hi')
# run_parallel_task(possible_states)


# # if __name__ == "__main__":
# #     with ProcessPoolExecutor() as executor:
# #         # Run 1000 tasks in parallel without passing any arguments
# #         futures = [executor.submit(task) for _ in range(500)]

# #     # Combine all returned sets into one using set union
# #     combined_states = set()

# #     # Iterate over the results from each future
# #     for future in futures:
# #         combined_states = combined_states.union(future.result())

# #     print(len(combined_states))

# #     with open('3pmap.pkl', 'wb') as f:
# #         pickle.dump(list(combined_states), f)


# for i in range(100):

#     mat= T*from_orbit([H,S,T]) * operator(5,1,[[e0],[n],[n],[n],[n]])

#     if mat.sde_sum()%5 == 0:
#         try:
#             print('got one')
#             print(mat.sde_profile())

#             m, o = synth_search(mat, full_set)
#             print(m.sde_profile())
#         except Exception:
#             print(mat)
#             print(mat.sde_profile())



#If we know the exact set of residues that clifford + R has then we are perfect. I could potentially upper bound this quantity and then try to collect them all, but thats iffy.

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