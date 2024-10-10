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

power_map = [e0, e1, e2]
def D_gate(a,b,c):
    return(operator(3,3, [[power_map[a],n,n], [n,power_map[b],n], [n,n,power_map[c]]]  ))

H.string = 'H'

R = operator(3,3,[[e0,n,n],[n,e0,n],[n,n,(-1)*e0]])

S = operator(3,3,[[e0,n,n],[n,e1,n],[n,n,e0]])
R.string = 'R'

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
    
    a= random.randint(0,2)
    if a == 0:
        return(go_stupid(argument*H, count+1, depth))
    elif a==1:
        return(go_stupid(argument*S, count+1, depth))
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
R_options = ['1','R']
D_options = ['1', 'D_gate(1,0,0)', 'D_gate(0,1,0)','D_gate(0,0,1)','D_gate(1,1,0)','D_gate(1,0,1)','D_gate(0,1,1)', 'D_gate(2,1,0)', 'D_gate(1,2,0)', 'D_gate(1,0,2)', 'D_gate(0,1,2)', 'D_gate(2,0,1)',' D_gate(0,2,1)']


all_options = [a + '*' + f for a in ['A', 'B', 'C', 'D'] for f in D_options]


def synth_search(oper):

    old_mat = oper


    for option in all_options:
        new_mat = eval(option) * old_mat
        if new_mat.sde < old_mat.sde:
            return(new_mat, option)
        
    
        
def column_sum(operator):
    matrix = operator.matrix
    sde = operator.sde
    sum_first_column = np.sum(matrix[:, 0])

    return(sum_first_column.sde+1< sde, sum_first_column.sde<sde)

def neighbors_mat(mat):
    A.string = 'H'
    B.string = 'HR'
    C.string = 'HRHH'
    D.string = 'HRHHR'
    neighbors = [A*mat, B*mat, C*mat, D*mat]

    return(neighbors)

AT = H
BT = R*H
CT = H*H*R*H
DT = R*H*H*R*H







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

mat= go_stupid()

print(mat)

while mat.sde > 1:
    mat, circ = synth_search(mat)
    print(circ)
    print(mat)

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


    



