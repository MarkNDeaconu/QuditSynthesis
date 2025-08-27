from quditsynthesis.datastructures import *
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

ket0 = operator(5,1, [
    [e0],
    [n],
    [n],
    [n],
    [n]
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

print(T)
total_list = [H,S,T,R]

# with open('cliffords5.pkl', 'rb') as f:
#     cliffords = pickle.load(f)

# for cliff in cliffords:
#     cliff.string = ''

# with open('cliffords5.pkl', 'wb') as f:
#     pickle.dump(cliffords, f)


# full_set = [h* D_gate(0,b,c,d,e) for b in range(5) for c in range(5) for d in range(5) for e in range(5) for h in z5.quotient(cliffords, z5.torus(cliffords, n))]

# cliffords = z5.subgroup([H,S])

'''cliffords = cliffords + [(-1)* c for c in cliffords]
print(len(cliffords))
print(len(z5.torus(cliffords, n)))

semi_cliff = z5.quotient(cliffords, z5.torus(cliffords, n))

print(len(semi_cliff))'''

# semi_cliff = z5.quotient(cliffords, z5.torus(cliffords, n))

# half_set = [H* D_gate(0,b,c,d,e) for b in range(5) for c in range(5) for d in range(5) for e in range(5)]

# print(len(half_set))

# r_modif_og = [ operator(5,5,[[ a*e0,n,n,n,n],[n,b*e0,n,n,n],[n,n,c*e0,n,n],[n,n,n,d*e0,n],[n,n,n,n,e*e0]]) for a in bin for b in bin for c in bin for d in bin for e in bin]
# r_modif = []
# for gate in r_modif_og:
#     if not((-1)*gate in r_modif):
#         r_modif.append(gate)

# print(len(r_modif))

# full_set = [a*b *c for a in half_set for b in half_set for c in r_modif]


# print(len(full_set))
# r_modif = [ operator(5,5,[[ e0,n,n,n,n],[n,b*e0,n,n,n],[n,n,c*e0,n,n],[n,n,n,d*e0,n],[n,n,n,n,e*e0]]) for b in bin for c in bin for d in bin for e in bin]


D_gates = [ D_gate(0,b,c,d,e) for b in range(5) for c in range(5) for d in range(5) for e in range(5)]
# full_set = [H* D_gate(0,b,c,d,e)* r for b in range(5) for c in range(5) for d in range(5) for e in range(5) for r in r_modif]
# full_set = [c* s for c in semi_cliff for s in D_gates]
# print(len(full_set))
#full_set = H*D*cliff*D
# full_set = [a * b  for a in [I,H*T,H*T*H*H*T] for b in cliffords]

# sde1set = []
# for elem in full_set:
#     if elem.sde ==1:
#         sde1set.append(elem)
# print(len(sde1set))

# H_options = [I, H, H*H , H*H*H]
# T_options = [I, T , T*T, T*T*T , T*T*T*T]

# all_options = [H*d*c for d in T_options for c in cliffords]

# print(len(all_options))
# options = [op for op in all_options if op.sde==1]

# mat = z5.from_orbit([H,S,T],200)

# print(mat.synthesize(options))


edges = [I,B,D]


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

        
# reduced_cyclotomics = []

# for v in range(5):
#     for w in  range(5):
#         for x in range(5):
#             for y in range(5):

#                 cyc = cyclotomic_element(z5, [v, w, x, y, 0],10)
#                 if cyc.sde == 10 and (cyc*cyc).sde == 20:
#                     reduced_cyclotomics.append(cyc)

# print(len(reduced_cyclotomics))

# cyclotomics_uptod = []
# for cyc in reduced_cyclotomics:
#     included = False
#     variants = [ (s * cyc * e1.power(k)).pmap_elem() for s in [-1,1] for k in range(5)]
#     for var in variants:
#         if var in cyclotomics_uptod:
#             included = True

#     if not(included):
#         cyclotomics_uptod.append(cyc)






# for i in range(100):
#     if (z5.from_orbit([H,S,R]) * ket0).synth_search(full_set) == None:
#         print('fail R')

# for i in range(10):
#     try:
#         mat = z5.from_orbit([H,S,R])
#         mat.synthesize(full_set)
#         print('yo')
#     except Exception:  
#         print('Fail synth R')
# for i in range(100):
#     if (z5.from_orbit([H,S,T]) * ket0).synth_search(full_set) == None:
#         print('fail')

# for i in range(10):
#     try:
#         mat = z5.from_orbit([H,S,T])
#         mat.synthesize(full_set)
#         print('yo')
#     except Exception:
#         print('Fail synth T')

# for i in range(100):
#     if (z5.from_orbit([H,S,T] + r_modif) * ket0).synth_search(full_set) == None:
#         print('fail full')

# for i in range(10):
#     try:
#         mat = z5.from_orbit([H,S,T] + r_modif)
#         mat.synthesize(full_set)
#         print('yo')
#     except Exception:
#         print('Fail full synth')
# for i in range(1000):
#     if (z5.from_orbit([H,S,R]) * ket0).synth_search(full_set) == None:
#         print('fail R')

# print(len(z5.quotient(cliffords, z5.torus(cliffords, n))))



# H.string = ''
# T.string = ''
# S.string = ''
# R.string = ''
# for e in cyclotomics_uptod:
#     possible_state = operator(5,1,[[e0],[e0],[e0],[e0],[e]]).rand_search([H,S,T,R])
#     print('done')

# set1 = [I, H*R, H*R*H*H*R]

# full_set = [s1 *s2 for s1 in set1 for s2 in cliffords]
# print(len(full_set))
# print(len(set(full_set)))
# print(z5.from_orbit([H,S,T,R]).sde_sum())
# for i in range(10):
#     print(z5.from_orbit([H,S,T,R]).synth_search(full_set))
# for sta in possible_states:
#     if sta.synth_search(full_set) == None:
#         print(sta)

# print('done')



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
