from datastructures import *
import random
import pickle


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

def go_stupid(argument=H, count=0, depth = random.randint(170,200)):
    if count >depth:
        return(argument)
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth))
    else:
        return(go_stupid(argument*T, count+1, depth))


def synth_search(oper):
    H_options = ['1','H','H*H','H*H*H']
    T_options = ['1','T','T*T','T*T*T','T*T*T*T']

    old_mat = oper

    all_options = [e+'*'+d+'*'+c+'*'+b  for e in H_options for d in T_options for c in H_options for b in T_options]

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) < np.sum(old_mat.sde_profile()):
            return(new_mat, option)
        
    all_options = [e+'*'+d+'*'+c+'*'+b + '*' + a + '*'+ f  for e in H_options for d in T_options for c in H_options for b in T_options for a in H_options[1:] for f in T_options[1:]]

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) < np.sum(old_mat.sde_profile()):
            return(new_mat, option)
        
    
        
def column_sum(operator):
    matrix = operator.matrix
    sde = operator.sde
    sum_first_column = np.sum(matrix[:, 0])

    return(sum_first_column.sde+1< sde, sum_first_column.sde<sde)

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

with open('5ditmat.pkl', 'rb') as file:
    mat = pickle.load(file)
    mat0=mat

# with open('5ditmat.pkl', 'wb') as file:
#     mat = go_stupid()
#     pickle.dump(mat, file)


print(mat)
string = ''
while mat.sde >2:
    mat, new_string = synth_search(mat)
    print(mat.sde_profile())

    print('')

    print(new_string)

    print('')

    string = new_string+'*'+string


print(mat)

print(string)


# print(string.count('T'))

# print(mat.sde_profile())
# print('')
# print((H*mat).sde_profile())
# print('')
# print( (H*H*mat).sde_profile())
# print('')
# print((H*H*H*mat).sde_profile())
# print(column_sum(mat))

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