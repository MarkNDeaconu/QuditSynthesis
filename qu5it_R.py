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

R = operator(5,5,[
    [e0, n, n, n, n],  
    [n, e0, n, n, n], 
    [n, n, e0, n, n], 
    [n, n, n, e0, n], 
    [n, n, n, n, (-1)*e0],  
])



X = operator(5,5,[
    [n, n, n, n, e0],
    [e0, n, n, n, n],  
    [n, e0, n, n, n], 
    [n, n, e0, n, n], 
    [n, n, n, e0, n] ]) 

#mat = H*T*H*T*H*H*T*T*H*T*T*H*T*H*T*T*T*H*H*H*T
# mat= H*T*T*H*T*H*T*H*T*T*H*H*T*H*T*H*T*T*H*T*H*T*T*T*H*T*T*T*H*T*T

# print(mat)
# random.seed(10)

def go_stupid(argument=H, count=0, depth = random.randint(170,200), H_count = 1):
    if count >depth:
        return(argument)
    
    a= random.randint(0,1)
    if a:
        return(go_stupid(argument*H, count+1, depth, H_count+1))
    else:
        return(go_stupid(argument*R, count+1, depth, H_count))


H_options = ['1','H','H*H','H*H*H']
T_options = ['1','R']

all_options = [c+'*'+b + '*' + a + '*'+ f   for c in H_options for b in T_options for a in H_options for f in T_options]

def synth_search(oper):

    old_mat = oper

    for option in all_options:
        new_mat = eval(option) * old_mat
        if np.sum(new_mat.sde_profile()) < np.sum(old_mat.sde_profile()):
            return(new_mat, option)
        
    
        
def column_sum(operator):
    matrix = operator.matrix
    sde = operator.sde
    sum_first_column = np.sum(matrix[:, 0])

    return(sum_first_column.sde+1< sde, sum_first_column.sde<sde)


# with open('5ditmat.pkl', 'rb') as file:
#     mat = pickle.load(file)
#     mat0=mat


with open('5ditmat.pkl', 'wb') as file:
    mat = go_stupid()
    pickle.dump(mat, file)


# mat1 = mat*H
# mat2 = mat * R*H
# print([mat1.sde, mat.sde, mat2.sde])


# mat3 = mat*H*H*H
# mat4 = mat*R*H*H*H

# print([mat3.sde, mat.sde, mat4.sde])

# print((mat*H*H).sde)


# sde_change = set()
# for i in range(10):
#     mat = go_stupid()

#     sde_list = [(mat*H).sde, (mat).sde, (mat*R*H).sde]

#     compare_list = [(mat*H*H*H).sde, (mat).sde, (mat*R*H*H*H).sde]

#     if sde_list != compare_list:
#         print('hi')


#     sde_change.add(((mat*H).sde - min(sde_list), (mat).sde - min(sde_list), (mat*R*H).sde- min(sde_list)))

# print(sde_change)

# print(mat1)

# print(mat2)


dropping_set = []

# string = ''
# while mat.sde >2:
#     mat, new_string = synth_search(mat)
#     print(mat.sde_profile())

#     print('')

#     print(new_string)

#     dropping_set.add(new_string)

#     print('')

#     string = new_string+'*'+string

# print(string)

# print(dropping_set)
for i in range(100):
    mat = go_stupid()
    res = [mat.sde , (H*mat).sde, (H*R*mat).sde, (H*R*H*H*R*mat).sde, (H*R*H*H*mat).sde]

    new_res = (mat.sde - min(res) , (H*mat).sde - min(res) ,(H*R*mat).sde - min(res), (H*R*H*H*R*mat).sde - min(res), (H*R*H*H*mat).sde - min(res) )
    dropping_set.append(new_res)

print(dropping_set)

h_count = dropping_set.count((1,0,2,2,2))
hr_count = dropping_set.count((1,2,0,2,2))
hrhhr_count = dropping_set.count((1,2,2,0,2))
hrhh_count = dropping_set.count((1,2,2,2,0))

print(h_count)
print(hr_count)
print(hrhhr_count)
print(hrhh_count)
print('')
print(len(dropping_set) - h_count - hr_count - hrhh_count - hrhhr_count)
# print(H*R*H*H*R*H*R*H*H*R*H*R*H*H*R*H*R*H*H*1*H*R*H*H*1*H*R*H*H*R*1*1*H*R*H*R*H*H*R*H*R*H*H*1*1*1*H*R*1*1*H*R*H*R*H*H*1*H*R*H*H*1*H*R*H*H*1*1*1*H*1*mat0)()