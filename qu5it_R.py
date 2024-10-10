from datastructures import *
import random
import pickle
import itertools

z5 =  cyclotomic_ring(5,math.sqrt(5))

n = cyclotomic_element(z5, [0,0,0,0,0]) 
e0 = cyclotomic_element(z5, [1,0,0,0,0])
e1 = cyclotomic_element(z5, [0,1,0,0,0])
e2 = cyclotomic_element(z5, [0,0,1,0,0])
e3 = cyclotomic_element(z5, [0,0,0,1,0])
e4 = cyclotomic_element(z5, [0,0,0,0,1])

power_map = [e0, e1, e2, e3, e4]

def D_gate(a,b,c,d,e):
    return(operator(5,5, [[power_map[a],n,n,n,n], [n,power_map[b],n,n,n], [n,n,power_map[c],n,n], [n,n,n,power_map[d],n], [n, n, n, n, power_map[e]]]))




H = (1/math.sqrt(5))*operator(5,5, [[e0,e0,e0,e0,e0], [e0,e1,e2,e3,e4], [e0,e2,e4,e1,e3], [e0,e3,e1,e4,e2], [e0, e4, e3,e2,e1]])
H.string = 'H'
R = operator(5,5,[
    [e0, n, n, n, n],  
    [n, e0, n, n, n], 
    [n, n, e0, n, n], 
    [n, n, n, e0, n], 
    [n, n, n, n, (-1)*e0],  
])


R.string = 'R'

T = operator(5,5, [[e0,n,n,n,n],[n,e1,n,n,n],[n,n,e3,n,n],[n,n,n,e2,n],[n,n,n,n,e4]])

S= operator(5,5, [[e0,n,n,n,n],[n,e1,n,n,n],[n,n,e3,n,n],[n,n,n,e1,n],[n,n,n,n,e0]])

print(S)
X = operator(5,5,[
    [n, n, n, n, e0],
    [e0, n, n, n, n],  
    [n, e0, n, n, n], 
    [n, n, e0, n, n], 
    [n, n, n, e0, n] ]) 

A= H
B= H*R
C = H*R*H*H
D= H*R*H*H*R
#mat = H*T*H*T*H*H*T*T*H*T*T*H*T*H*T*T*T*H*H*H*T
# mat= H*T*T*H*T*H*T*H*T*T*H*H*T*H*T*H*T*T*H*T*H*T*T*T*H*T*T*T*H*T*T

# print(mat)
# random.seed(10)


def go_stupid(argument=H, count=0, depth = random.randint(170,200), H_count = 1):
    if count >depth:
        return(argument)
    
    a= random.randint(0,2)
    if a == 1:
        return(go_stupid(argument*H, count+1, depth, H_count+1))
    elif a ==2:
        return(go_stupid(argument*R, count+1, depth, H_count))
    else:
        return(go_stupid(argument*S, count+1, depth, H_count))


H_options = ['1','H','H*H', 'H*H*H']
R_options = ['1','R']



D_options = ['D_gate(0,'+str(i) +', ' +str(j) + ',' + str(k) + ',' + str(l) + ')' for l in range(5) for i in range(5)for j in range(5)for j in range(5) for k in range(5)]



all_options = [a + '*' + f for a in ['A', 'B', 'C', 'D'] for f in D_options]


def synth_search(oper):

    old_mat = oper


    for option in all_options:
        new_mat = eval(option) * old_mat
        if new_mat.sde < old_mat.sde:
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


def make_hashable(matrix):
    # Recursively convert all lists to tuples
    if isinstance(matrix, list):
        return tuple(make_hashable(element) for element in matrix)
    return matrix


def neighbors_mat(mat):
    A.string = 'H'
    B.string = 'HR'
    C.string = 'HRHH'
    D.string = 'HRHHR'
    neighbors = [A*mat, B*mat, C*mat, D*mat]
    return(neighbors)

# listy = set()

# for i in range(200):
#     pmap_set = [[set(),set(),set(), set(), set()],[set(),set(),set(), set(), set()],[set(),set(),set(), set(), set()],[set(),set(),set(), set(), set()],[set(),set(),set(), set(), set()]]
#     mat = go_stupid()
#     pmaps = mat.pmap()
#     for a in range(5):
#         for b in range(5):
#             pmap_set[a][b] = make_hashable(pmaps[a][b])
#     listy.add(make_hashable(pmap_set))

# print(len(listy))

# print(len(pmap_set[0][0]))
# print(len(pmap_set[0][1]))
# print(len(pmap_set[0][2]))
# print(len(pmap_set[0][3]))
# print(len(pmap_set[0][4]))
# print(len(pmap_set[1][0]))
# print(len(pmap_set[1][1]))
# print(len(pmap_set[1][2]))
# print(len(pmap_set[1][3]))
# print(len(pmap_set[1][4]))
# print(len(pmap_set[2][0]))
# print(len(pmap_set[2][1]))
# print(len(pmap_set[2][2]))
# print(len(pmap_set[2][3]))
# print(len(pmap_set[2][4]))
# print(len(pmap_set[3][0]))
# print(len(pmap_set[3][1]))
# print(len(pmap_set[3][2]))
# print(len(pmap_set[3][3]))
# print(len(pmap_set[3][4]))
# print(len(pmap_set[4][0]))
# print(len(pmap_set[4][1]))
# print(len(pmap_set[4][2]))
# print(len(pmap_set[4][3]))
# print(len(pmap_set[4][4]))



# five_list = [0,1,2,3,4]
# cart_prod = itertools.product(five_list,five_list,five_list,five_list,[0])
# mod_5_list = []
# for element in cart_prod:
#     mod_5_list.append(element)

# mod_5 = []
# for element in mod_5_list:
#     if z5.reduced(list(element), 0) ==(list(element),0):
#         mod_5.append(element)

# print(len(mod_5))

# print(z5.loc_char)
# mat = H*H*H*H
# print(mat*H)
# print(mat*H*R)

# print(mat*H*R*H*H)
# print(mat*H*R*H*H*R)



'''dropping_set = []

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
'''

mat = go_stupid()
# print(synth_search(mat))
# print(mat)
print(mat)
while mat.sde >1:
    mat, string = synth_search(mat)
    
    print(string)
    print(mat)


'''lowers = [A,B,C,D]
def can_it_be_reduced(mat):
    return(mat.sde-1 > min( [((x*mat).sde) for x in lowers for y in lowers]))'''

# print((A*mat).sde)
# print((B*mat).sde)
# print((C*mat).sde)
# print((D*mat).sde)

# print('')

# mat = S*mat
# print(mat.sde)
# print((A*mat).sde)
# print((B*mat).sde)
# print((C*mat).sde)
# print((D*mat).sde)



'''def find_candidate(mat):
    for i in range(4):
        for j in range(4):
            for m in range(4):
                for l in range(4):
                    for k in range(4):
                        if can_it_be_reduced(D_gate(i,j,m,l,k) * mat):
                            return(D_gate(i,j,m,l,k) * mat)
                        
print(find_candidate(mat))
'''

