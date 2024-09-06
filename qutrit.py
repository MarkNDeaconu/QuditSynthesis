from datastructures import *
import random
import pickle

from complex_verification import verify

z3 =  cyclotomic_ring(3,complex(0, math.sqrt(3)))
n = cyclotomic_element(z3, [0,0,0]) 
e0 = cyclotomic_element(z3, [1,0,0])
e1 = cyclotomic_element(z3, [0,1,0])
e2 = cyclotomic_element(z3, [0,0,1])


H = (1/math.sqrt(3))*operator(3,3, [[e0,e0,e0], [e0,e1,e2], [e0,e2,e1]])
I = H*H*H*H
R = operator(3,3,[[e0,n,n],[n,e0,n],[n,n,(-1)*e0]])

# print(H)
# print(H*H*H*H)

# print(H*H)

# print(z3.loc_char)

# print( cyclotomic_element(z3, [9,0,0]))

def go_stupid(argument=R, count=0, depth = random.randint(100,200)):
    if count >depth:
        return(argument)
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth))
    
    else:
        return(go_stupid(argument*R, count+1, depth))
    

def synth_search(oper):
    H_options = ['1','H','H*H','H*H*H']
    R_options = ['1','R']

    old_mat = oper

    all_options = [e+'*'+d+'*'+c+'*'+b  for e in H_options for d in R_options for c in H_options for b in R_options]

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) < np.sum(old_mat.sde_profile()):
            return(new_mat, option)

        
def column_sum(operator):
    matrix = operator.matrix
    sde = operator.sde
    sum_first_column = np.sum(matrix[:, 0])

    return(sum_first_column.sde+1< sde, sum_first_column.sde<sde)

# print(go_stupid(H))


print(H)
print(R)


with open('3ditmat.pkl', 'rb') as file:
    mat = pickle.load(file)


mat0=R*H*R*H*R*H*(H*H*H*R)

print(mat0)
print((mat0).sde_profile())


print(mat)

print(H.comp())
print(R.comp()*H.comp())



# with open('3ditmat.pkl', 'wb') as file:
#     mat = go_stupid()
#     pickle.dump(mat, file)

'''string = ''
while mat.sde >2:
    mat, new_string = synth_search(mat)
    print(mat.sde_profile())

    print('')

    print(new_string)

    print('')

    string = new_string+'*'+string


print(mat)

print(string)'''
