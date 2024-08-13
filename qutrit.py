from datastructures import *
import random
z3 =  cyclotomic_ring(3,math.sqrt(3))
n = cyclotomic_element(z3, [0,0,0]) 
e0 = cyclotomic_element(z3, [1,0,0])
e1 = cyclotomic_element(z3, [0,1,0])
e2 = cyclotomic_element(z3, [0,0,1])



H = (1/math.sqrt(3))*operator(3,3, [[e0,e0,e0], [e0,e1,e2], [e0,e2,e1]])

T = operator(3,3,[[e1,n,n],[n,e0,n],[n,n,e2]])

# print(H)
# print(H*H*H*H)

# print(H*H)

# print(z3.loc_char)

# print( cyclotomic_element(z3, [9,0,0]))

def go_stupid(argument, count=0, depth = random.randint(100,200)):
    if count >depth:
        return(argument)
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth))
    
    else:
        return(go_stupid(argument*T, count+1, depth))
    

# print(go_stupid(H))

print(H*T*H*T*T*T*H*H*T*H*T*T*T)
print(T*T*H*T*H*T*T*T*H*H*T*H*T*T*T*T)