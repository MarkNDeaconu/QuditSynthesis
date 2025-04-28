from datastructures import *
import random
import pickle
import itertools
from concurrent.futures import ProcessPoolExecutor


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

with open('QuditSynthesis/cliffords3.pkl', 'rb') as f:
    cliffords = pickle.load(f)

print(len(cliffords))

print(len(z3.torus(cliffords, n)))

print(D.monomial_check())


full_set = [a * b  for a in [I,B] for b in cliffords]

edges = [A,B,C,D]


orbit = z3.subgroup_bfs([H,S,R], 8)

mono = set()
for elem in orbit:
    if elem.monomial_check():
        mono.add(elem)

print(len(mono))

# def task():
#     curr = H
#     states = set()
#     for i in range(3000):
#         if curr.sde > 3:
#             states.update(set(curr.pmap_state()))
#         curr = random.choice([H,S,R]) * curr

#     return(states)


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
    


# # with open('3pmap.pkl', 'rb') as f:
# #     results = pickle.load(f)

# all_cyclotomics = [cyclotomic_element(z3, [x, y, 0],10) for x in range(3) for y in range(3)]
# reduced_cyclotomics = []

# for cyc in all_cyclotomics:
#     if cyc.sde == 10:
#         reduced_cyclotomics.append(cyc)

# all_states = [operator(3,1,[[a],[b],[c]]) for a in reduced_cyclotomics for b in reduced_cyclotomics for c in reduced_cyclotomics]

# print(len(all_states))


# # new_set = []

# # for x in results:
# #     if x.sde ==10:
# #         new_set.append(x)

# for res in all_states:
#     try:
#         synth_search(res, full_set)
#     except Exception:
#         print('fail')

