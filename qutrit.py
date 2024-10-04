from datastructures import *
import random
import pickle
import itertools
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

R = operator(3,3,[[e0,n,n],[n,e0,n],[n,n,(-1)*e0]])

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R
# print(H)
# print(H*H*H*H)

# print(H*H)

# print(z3.loc_char)

# print( cyclotomic_element(z3, [9,0,0]))

def go_stupid(argument=A, count=0, depth = random.randint(130,140), string = ''):
    if count >depth:
        return(argument)
    
    a= random.randint(0,1)
    if a == 0:
        return(go_stupid(argument*H, count+1, depth))
    else:
        return(go_stupid(argument*R, count+1, depth))
    # a= random.randint(0,3)
    # if a == 1:
    #     return(go_stupid(argument*A, count+1, depth, string + 'A'))
    
    # elif a==2 :
    #     return(go_stupid(argument*B, count+1, depth, string + 'B'))
    
    # elif a==3:
    #     return(go_stupid(argument*C, count+1, depth, string + 'C'))
    # else:
    #     return(go_stupid(argument*D, count+1, depth, string + 'D'))
    

H_options = ['1','H','H*H','H*H*H']
T_options = ['1','R']

#all_options = [c+'*'+b + '*' + a + '*'+ f   for c in H_options for b in T_options for a in H_options for f in T_options]

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

def neighbors_mat(mat):
    neighbors = [(A*mat, 'H'), (B*mat, 'HR'), (C*mat, 'HRHH'), (D*mat,'HRHHR')]
    sorted_neigh = sorted(neighbors, key=lambda x: x[0])
    return(sorted_neigh)

AT = H
BT = R*H
CT = H*H*R*H
DT = R*H*H*R*H

mat= go_stupid()
print(mat)

for row in mat.pmap():
    print(row)

# print((A*mat*AT).sde)
# print((A*mat*BT).sde)
# print((A*mat*CT).sde)
# print((A*mat*DT).sde)
# print((B*mat*AT).sde)
# print((B*mat*BT).sde)
# print((B*mat*CT).sde)
# print((B*mat*DT).sde)
# print((D*mat*AT).sde)
# print((D*mat*BT).sde)
# print((D*mat*CT).sde)
# print((D*mat*DT).sde)
# print((C*mat*AT).sde)
# print((C*mat*BT).sde)
# print((C*mat*CT).sde)
# print((C*mat*DT).sde)





'''
with open('5ditmat.pkl', 'wb') as file:
    mat = go_stupid()
    pickle.dump(mat, file)

print(mat)


print(mat.sde_profile())




string = ''
while mat.sde >2:
    mat, new_string = synth_search(mat)
    print(mat.sde_profile())

    print('')

    print(new_string)

    print('')

    string = new_string+'*'+string

print(string)
'''
print(mat)



print(H*A)
print(H*B)
print(H*C)
print(H*D)