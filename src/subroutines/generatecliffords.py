import random
import pickle
from quditsynthesis.datastructures import *

def generate_cliffords(ring: cyclotomic_ring, dim, operator_list, depth = 10):
    
    cliffords = ring.subgroup_bfs(operator_list, depth)

    print(len(cliffords))
    with open('stored_sets.cliffords'+str(dim) + '.pkl', 'wb') as f:
        pickle.dump(cliffords, f)
    
    


# with open('cliffords3.pkl', 'rb') as f:
#     cliffords = pickle.load(f)
#     cliffords = list(cliffords)


# final_cliffords = []

# for c in cliffords:
#     if not((-1) * c in final_cliffords):
#         c.string = ''
#         final_cliffords.append(c)


# with open('cliffords3.pkl', 'wb') as f:
#     pickle.dump(final_cliffords, f)

# print(len(final_cliffords))

