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

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R


def go_stupid(argument=R, count=0, depth = random.randint(100,200)):
    if count >depth:
        return(argument)
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth))
    
    else:
        return(go_stupid(argument*R, count+1, depth))


H_options = ['1','H','H*H','H*H*H']
R_options = ['1','R']

all_options = [e+'*'+d+'*'+c+'*'+b  for e in H_options for d in R_options for c in H_options for b in R_options]

def synth_search(oper):

    old_mat = oper

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) < np.sum(old_mat.sde_profile()):
            return(new_mat, option)
        
def neighbors_mat(mat):
    A.string = 'H'
    B.string = 'HR'
    C.string = 'HRHH'
    D.string = 'HRHHR'
    neighbors = [A*mat, B*mat, C*mat, D*mat]
    return(neighbors)


"""mat= go_stupid()

print(mat.sde_profile())

print(mat)

string = ''
while mat.sde >2:
    mat, new_string = synth_search(mat)
    print(mat.sde_profile())

    print('')

    print(new_string)

    print('')

    string = new_string+'*'+string

print(string)"""

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

print(H*H)