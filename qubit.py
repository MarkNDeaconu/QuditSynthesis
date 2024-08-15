from datastructures import *
import random
z8 =  cyclotomic_ring(8,math.sqrt(2))
n = cyclotomic_element(z8, [0,0,0,0,0,0,0,0]) 
e0 = cyclotomic_element(z8, [1,0,0,0,0,0,0,0])
e1 = cyclotomic_element(z8, [0,1,0,0,0,0,0,0])
e2 = cyclotomic_element(z8, [0,0,1,0,0,0,0,0])
e3 = cyclotomic_element(z8, [0,0,0,1,0,0,0,0])

H = (1/math.sqrt(2))*operator(2,2, [[e0,e0], [e0,(-1)*e0]])

T= operator(2,2,[[e0,n],[n, e1]])

print(H*H*H*T)

print(T.power(8))

def random_mat(operator= H, depth = 50):
    new_oper = operator
    for op in range(depth):
        a = random.randint(1,8)
        new_oper = H*T.power(a) *new_oper

    return(new_oper)

# def synthesize(operator, target = 3):



print(random_mat())