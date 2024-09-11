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



T = operator(7,7,[
    [e0, n, n, n, n, n, n],  
    [n, e0, n, n, n, n, n], 
    [n, n, e0, n, n, n, n], 
    [n, n, n, e0, n, n, n], 
    [n, n, n, n, e0, n, n],  
    [n, n, n, n, n, e0, n],  
    [n, n, n, n, n, n, (-1)*e0]  
])

# T2 = T*T
# T3 = T*T*T
# T4 = T*T*T*T
# T5 = T*T*T*T*T
# T6 = T*T*T*T*T*T

def go_stupid(argument=T, count=0, depth = random.randint(100,200)):
    if count >depth:
        return(argument)
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth))
    
    else:
        return(go_stupid(argument*T, count+1, depth))

def raise_sde(oper):
    H_options = ['1','H','H2','H3']
    T_options = ['1','T']

    old_mat = oper

    all_options = [e+'*'+d+'*'+c+'*'+b + '*' +f + '*' + g   for e in H_options for d in T_options for c in H_options for b in T_options for f in H_options for g in T_options]

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) > np.sum(old_mat.sde_profile()):
            return(new_mat, option)
    

def synth_search(oper):
    H_options = ['1','H','H*H','H*H*H']
    T_options = ['1','T','T*T','T*T*T','T*T*T*T', 'T*T*T*T*T', 'T*T*T*T*T*T']

    old_mat = oper

    all_options = [e+'*'+d+'*'+c+'*'+b  for e in H_options for d in T_options for c in H_options for b in T_options]

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) < np.sum(old_mat.sde_profile()):
            return(new_mat, option)
        

mat= go_stupid()

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

print(string)
