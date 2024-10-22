from datastructures import *
import random
import pickle


z7 =  cyclotomic_ring(7,complex(0, math.sqrt(7)))

n = cyclotomic_element(z7, [0,0,0,0,0,0,0]) 
e0 = cyclotomic_element(z7, [1,0,0,0,0,0,0])
e1 = cyclotomic_element(z7, [0,1,0,0,0,0,0])
e2 = cyclotomic_element(z7, [0,0,1,0,0,0,0])
e3 = cyclotomic_element(z7, [0,0,0,1,0,0,0])
e4 = cyclotomic_element(z7, [0,0,0,0,1,0,0])
e5 = cyclotomic_element(z7, [0,0,0,0,0,1,0])
e6 = cyclotomic_element(z7, [0,0,0,0,0,0,1])

e0m = cyclotomic_element(z7, [1,0,0,0,0,0,0],1)
e1m = cyclotomic_element(z7, [0,1,0,0,0,0,0],1)
e2m = cyclotomic_element(z7, [0,0,1,0,0,0,0],1)
e3m = cyclotomic_element(z7, [0,0,0,1,0,0,0],1)
e4m = cyclotomic_element(z7, [0,0,0,0,1,0,0],1)
e5m = cyclotomic_element(z7, [0,0,0,0,0,1,0],1)
e6m = cyclotomic_element(z7, [0,0,0,0,0,0,1],1)





H = operator(7,7, [
    [e0m, e0m, e0m, e0m, e0m, e0m, e0m],      
    [e0m, e1m, e2m, e3m, e4m, e5m, e6m],     
    [e0m, e2m, e4m, e6m, e1m, e3m, e5m],     
    [e0m, e3m, e6m, e2m, e5m, e1m, e4m],      
    [e0m, e4m, e1m, e5m, e2m, e6m, e3m],     
    [e0m, e5m, e3m, e1m, e6m, e4m, e2m],     
    [e0m, e6m, e5m, e4m, e3m, e2m, e1m]       
])

H2 = H*H
H3 = H*H*H



R = operator(7,7,[
    [e0, n, n, n, n, n, n],  
    [n, e0, n, n, n, n, n], 
    [n, n, e0, n, n, n, n], 
    [n, n, n, e0, n, n, n], 
    [n, n, n, n, e0, n, n],  
    [n, n, n, n, n, e0, n],  
    [n, n, n, n, n, n, (-1)*e0]  
])

I= R*R

S = operator(7,7,[
    [e0, n, n, n, n, n, n],  
    [n, e1, n, n, n, n, n], 
    [n, n, e3, n, n, n, n], 
    [n, n, n, e6, n, n, n], 
    [n, n, n, n, e3, n, n],  
    [n, n, n, n, n, e1, n],  
    [n, n, n, n, n, n, e0]  
])

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R

B.string = 'B'
D.string = 'D'


def go_stupid(argument=H, count=0, depth = random.randint(100,200), H_count = 1):
    if count >depth:
        return(argument)
    
    a= random.randint(0,2)
    if a == 1:
        return(go_stupid(argument*H, count+1, depth, H_count+1))
    elif a ==2:
        return(go_stupid(argument*R, count+1, depth, H_count))
    else:
        return(go_stupid(argument*S, count+1, depth, H_count))


H_options = ['1','H','H*H','H*H*H']
R_options = ['1','R']

# with open('cliffords7.pkl', 'rb') as f:
#     cliffords = pickle.load(f)

# full_set = [a * b  for b in cliffords for a in [I,B,D]]

with open('cliffords7fullset.pkl', 'rb') as f:
    full_set = pickle.load(f)

# with open('cliffords7fullset.pkl', 'wb') as file:
#     pickle.dump(full_set, file)

# semi_set = [a * c * d  for d in H_options for c in S_options for a in dropping]

def synth_search(oper):

    old_mat = oper

    for option in full_set:
        new_mat = option * old_mat
        if new_mat.sde < old_mat.sde:
            return(new_mat, option.string)


def neighbors_mat(mat):
    A.string = 'H'
    B.string = 'HR'
    C.string = 'HRHH'
    D.string = 'HRHHR'
    neighbors = [A*mat, B*mat, C*mat, D*mat]
    return(neighbors)


for i in range(100):
    mat = go_stupid()



    while mat.sde >1:
        mat, string = synth_search(mat)
        print(mat)
        print(string)
    print(i)


# cliffords = set()
# curr_mat = H

# for i in range(5000000):
#     cliffords.add(curr_mat)

#     a= random.randint(0,3)

#     if a ==0:
#         curr_mat = H*curr_mat
#     elif a== 1:
#         curr_mat = S*curr_mat
#     elif a==2:
#         curr_mat = H*H*curr_mat
#     else:
#         curr_mat = S*S*curr_mat

# print(len(cliffords))

# with open('cliffords7.pkl', 'wb') as file:
#     pickle.dump(list(cliffords), file)

