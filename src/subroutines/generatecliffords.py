import random
import pickle
from quditsynthesis.datastructures import *
import os

import subprocess
ROOT = subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"]
).decode().strip()



def generate_cliffords(ring: cyclotomic_ring, dim, operator_list, depth = 10):
    
    cliffords = ring.subgroup_bfs(operator_list, depth)

    print(len(cliffords))

    with open(ROOT+'/src/stored_sets/cliffords'+str(dim) + '.pkl', 'wb') as f:
        pickle.dump(cliffords, f)

def store_group(group, group_name):

    with open(ROOT+'/src/stored_sets/'+group_name+'.pkl', 'wb') as f:
        pickle.dump(group, f)

def retrieve_group(group_name):

    with open(ROOT+'/src/stored_sets/'+group_name+'.pkl', 'rb') as f:
        group = pickle.load(f)
    
    return(group)

    
    


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

