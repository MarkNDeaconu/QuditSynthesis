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

def random_mat(operator= H, depth = 500):
    new_oper = operator
    for op in range(depth):
        a = random.randint(1,8)
        new_oper = H*T.power(a) *new_oper

    return(new_oper)

def synthesize(operator, target = 3):
    new_operator = operator
    sde = operator.sde
    string = ''
    while sde >target:
        if [(operator.matrix[0,0].ring.pmap(new_operator.matrix[0,0].coefficients)[i] + operator.matrix[0,0].ring.pmap(new_operator.matrix[1,0].coefficients)[i])%2 for i in range(4)] in [[0,0,0,0], [1,1,1,1]]:
            new_operator = T*H*new_operator
            string = 'T*H*' + string
        else:
            new_operator = T*new_operator
            string = 'T*' + string
        sde = new_operator.sde
    return(new_operator, string)


mat = random_mat()

print(mat)
print(synthesize(mat))
