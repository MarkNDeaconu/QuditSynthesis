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

power_map = [e0, e1, e2, e3, e4]

def D_gate(a,b,c,d,e):
    gate = operator(5,5, [[power_map[a],n,n,n,n], [n,power_map[b],n,n,n], [n,n,power_map[c],n,n], [n,n,n,power_map[d],n], [n, n, n, n, power_map[e]]])
    gate.string = 'D('+str(a)+str(b)+str(c)+str(d)+str(e)+')'
    return(gate)




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

T = operator(5,5, [[e0,n,n,n,n],[n,e1,n,n,n],[n,n,e3,n,n],[n,n,n,e2,n],[n,n,n,n,e4]])

S= operator(5,5, [[e0,n,n,n,n],[n,e1,n,n,n],[n,n,e3,n,n],[n,n,n,e1,n],[n,n,n,n,e0]])

S.string = 'S'

T.string = 'T'

X = operator(5,5,[
    [n, n, n, n, e0],
    [e0, n, n, n, n],  
    [n, e0, n, n, n], 
    [n, n, e0, n, n], 
    [n, n, n, e0, n] ]) 

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R


I = R*R
I.string = ''

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

def go_stupid2(argument=H, count=0, depth = random.randint(50,100), H_count = 1):
    if count >depth:
        return(argument)
    S.string = 'S'

    T.string = 'T'

    H.string = 'H'

    a= random.randint(0,2)
    if a == 1:
        return(go_stupid2(argument*H, count+1, depth, H_count+1))
    elif a ==2:
        return(go_stupid2(argument*T, count+1, depth, H_count))
    else:
        return(go_stupid2(argument*S, count+1, depth, H_count))


H_options = [I, H,H*H, H*H*H]
R_options = [I, R]

S_options = [I, S, S*S, S*S*S, S*S*S*S]

dropping = [A,B,C,D]





with open('cliffords5.pkl', 'rb') as f:
    cliffords = pickle.load(f)



full_set = [a * b  for a in [I,B,D] for b in cliffords]


# def synth_search(oper):
#     old_mat = oper
#     results = set()
#     for option in full_set:
#         new_mat = option * old_mat
#         if new_mat.sde < old_mat.sde:
#             results.add(option.string)
#     return(results)
I.string = 'C'
def synth_search(oper):
    old_mat = I*oper
    for option in full_set:
        new_mat = option * old_mat
        if new_mat.sde < old_mat.sde:
            return(new_mat, option.string)




        
    
        
def column_sum(operator):
    matrix = operator.matrix
    sde = operator.sde
    sum_first_column = np.sum(matrix[:, 0])

    return(sum_first_column.sde+1< sde, sum_first_column.sde<sde)


def neighbors_mat(mat):
    B.string = "HR"
    D.string = 'HRHHR'
    lowest_neighbor, option = synth_search(mat)
    if option == '':
        neighbors = [lowest_neighbor, B*mat, D*mat]
    elif option == 'HR':
        neighbors = [mat, lowest_neighbor,  D*mat]
    else:
         neighbors = [mat,  B*mat, lowest_neighbor]
    
    return(neighbors)



# for i in range(100):
#     mat0 = go_stupid2()
#     print(mat0.sde_profile())
#     mat = go_stupid() *mat0
#     print(mat.sde_profile())
#     try:
#         while mat.sde >1:
#             mat, string = synth_search(mat)
#             print(string)
#             print(mat.sde_profile())
#     except Exception:
#         print('fail')
#     print(i)



