from datastructures import *
import random
import pickle

from complex_verification import verify

z3 =  cyclotomic_ring(3,complex(0, math.sqrt(3)))
n = cyclotomic_element(z3, [0,0,0]) 
e0 = cyclotomic_element(z3, [1,0,0])
e1 = cyclotomic_element(z3, [0,1,0])
e2 = cyclotomic_element(z3, [0,0,1])

e0m = cyclotomic_element(z3, [1,0,0],1)
e1m = cyclotomic_element(z3, [0,1,0],1)
e2m = cyclotomic_element(z3, [0,0,1],1)

H = operator(3,3, [[e0m,e0m,e0m], [e0m,e1m,e2m], [e0m,e2m,e1m]])
I = H*H*H*H
R = operator(3,3,[[e0,n,n],[n,e0,n],[n,n,(-1)*e0]])



def go_stupid(argument=R, count=0, depth = random.randint(200,300)):
    if count >depth:
        return(argument)
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth))
    
    else:
        return(go_stupid(argument*R, count+1, depth))
    

H_options = ['1','H','H*H','H*H*H']
T_options = ['1','R']

all_options = [c+'*'+b + '*' + a + '*'+ f   for c in H_options for b in T_options for a in H_options for f in T_options]

def synth_search(oper):

    old_mat = oper

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) < np.sum(old_mat.sde_profile()):
            return(new_mat, option)
        
    
        
def column_sum(operator):
    matrix = operator.matrix
    sde = operator.sde
    sum_first_column = np.sum(matrix[:, 0])

    return(sum_first_column.sde+1< sde, sum_first_column.sde<sde)



with open('5ditmat.pkl', 'wb') as file:
    mat = go_stupid()
    pickle.dump(mat, file)

# with open('5ditmat.pkl', 'rb') as file:
#     mat = pickle.load(file)
#     mat0=mat

print(mat)

print(H*mat)

print(H*R*mat)

print(H*R*H*H*R*mat)

print(H*R*H*H*mat)

# dropping_set = []


# for i in range(100):
#     mat = go_stupid()
#     res = [mat.sde , (H*mat).sde, (H*R*mat).sde, (H*R*H*H*R*mat).sde, (H*R*H*H*mat).sde] 

#     new_res = (mat.sde - min(res) , (H*mat).sde - min(res) ,(H*R*mat).sde - min(res), (H*R*H*H*R*mat).sde - min(res), (H*R*H*H*mat).sde - min(res) )
#     dropping_set.append(new_res)

# print(dropping_set)

# h_count = dropping_set.count((1,0,2,2,2))
# hr_count = dropping_set.count((1,2,0,2,2))
# hrhhr_count = dropping_set.count((1,2,2,0,2))
# hrhh_count = dropping_set.count((1,2,2,2,0))

# print(h_count)
# print(hr_count)
# print(hrhhr_count)
# print(hrhh_count)
# print('')
# print(len(dropping_set) - h_count - hr_count - hrhh_count - hrhhr_count)


# print(H)

# print(H*R)

# print(H*R*H*H*R)
# print(H*R*H*H)

# print(z3.loc_char)