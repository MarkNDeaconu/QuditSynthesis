from datastructures import *
import random

z5 =  cyclotomic_ring(5,math.sqrt(5))

n = cyclotomic_element(z5, [0,0,0,0,0]) 
e0 = cyclotomic_element(z5, [1,0,0,0,0])
e1 = cyclotomic_element(z5, [0,1,0,0,0])
e2 = cyclotomic_element(z5, [0,0,1,0,0])
e3 = cyclotomic_element(z5, [0,0,0,1,0])
e4 = cyclotomic_element(z5, [0,0,0,0,1])



H = (1/math.sqrt(5))*operator(5,5, [[e0,e0,e0,e0,e0], [e0,e1,e2,e3,e4], [e0,e2,e4,e1,e3], [e0,e3,e1,e4,e2], [e0, e4, e3,e2,e1]])

T = operator(5,5, [[e0,n,n,n,n],[n,e1,n,n,n],[n,n,e3,n,n],[n,n,n,e2,n],[n,n,n,n,e4]])

I = H*H*H*H

#mat = H*T*H*T*H*H*T*T*H*T*T*H*T*H*T*T*T*H*H*H*T
# mat= H*T*T*H*T*H*T*H*T*T*H*H*T*H*T*H*T*T*H*T*H*T*T*T*H*T*T*T*H*T*T

# print(mat)
# random.seed(10)

def go_stupid(argument, count=0, depth = random.randint(100,200)):
    if count >depth:
        return(argument)
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth))
    else:
        return(go_stupid(argument*T, count+1, depth))


def synth_search(oper):
    H_options = [I,H,H*H,H*H*H]
    T_options = [I,T,T*T,T*T*T,T*T*T*T]

    all_options = [d*c*b*a for a in H_options for b in T_options for c in H_options for d in T_options]

    for option in all_options:
        new_mat = option * oper
        if new_mat.sde < oper.sde or new_mat.sde2 < oper.sde2:
            return(new_mat)

'''test = H*T*H*T*T*T*H*H*H*T*H*H*T*T*T*H*T*H*T*H*T*T*H*H*H*T
print(test)
print(T*T*T*T*H*H*H*test)
'''
'''test = go_stupid(H)

print(test)
final = synth_search(test)
print(final)

print(test.sde_profile())
print(final.sde_profile())'''

print(H.power(2))

# print(test)
# print(test.sde_profile())

# better = test

# print(better.sde_profile())
# print(better)

# print(z5.reduced(np.sum((better.matrix)[:, 0]).coefficients,0))


#print(test.comp()) 
#([-25, -8, -26, -1, 0], 0)

# print(H*H*H)

# print(T)
# print(T*T*T*T)