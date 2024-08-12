from datastructures import *

z8 =  cyclotomic_ring(8,math.sqrt(2))
n = cyclotomic_element(z8, [0,0,0,0,0,0,0,0]) 
e0 = cyclotomic_element(z8, [1,0,0,0,0,0,0,0])
e1 = cyclotomic_element(z8, [0,1,0,0,0,0,0,0])
e2 = cyclotomic_element(z8, [0,0,1,0,0,0,0,0])
e3 = cyclotomic_element(z8, [0,0,0,1,0,0,0,0])

H = (1/math.sqrt(2))*operator(2,2, [[e0,e0], [e0,(-1)*e0]])

T= operator(2,2,[[e0,n],[n, e1]])

print(H*H)


print(circulant(gauss_sequence(8)))