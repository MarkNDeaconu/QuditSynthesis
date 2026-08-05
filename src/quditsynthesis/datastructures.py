import math 
import numpy as np
from typing import Optional
from tabulate import tabulate
import random

superscript_map = {
    '0': '⁰',
    '1': '¹',
    '2': '²',
    '3': '³',
    '4': '⁴',
    '5': '⁵',
    '6': '⁶',
    '7': '⁷',
    '8': '⁸',
    '9': '⁹',
    '-': '⁻'
}

def circulant(row):
    n = len(row)
    circ_matrix = np.array([np.roll(row, i) for i in range(n)])
    return(circ_matrix)

def gauss_sequence(p):
    sequence = [0]*p
    for num in range(p):
        sequence[(num**2)%p] += 1

    return(sequence)

def multiply_many(operators):
    result = operators[0]
    for op in operators[1:]:
        result = result * op
    return(result)

def multiply_selected(generators, indices):
    acc = generators[indices[0]]
    out = [acc]
    for i in indices[1:]:
        acc = generators[i] * acc
        out.append(acc)
    return(out)


class cyclotomic_ring:
    def __init__(self, root_of_unity, localization) -> None:
        self.root_of_unity = root_of_unity
        self.localization = localization
        self.num_coefficient = root_of_unity
        # loc_char = p * circulant(gauss_sequence(p)): integer divisibility-test matrix for reduction
        self.loc_char = root_of_unity*circulant(gauss_sequence(root_of_unity))

            
    def __eq__(self, value: object) -> bool:
        if self.root_of_unity == value.root_of_unity and self.localization == value.localization:
            return True
        else:
            return False
    
    def add(self,coefficients1, coefficients2):
        new_val = []
        for i in range(self.num_coefficient):
            new_val.append(coefficients1[i] + coefficients2[i])
        return new_val

    def mul(self,coefficients1, coefficients2):
        new_val = [0] * self.num_coefficient
        for i in range(self.num_coefficient):
            for j in range(self.num_coefficient):
                new_val[(i+j)%self.num_coefficient]+=(coefficients1[i]*coefficients2[j])
        return new_val
    
    def matrix(self, coeff, matrix):
        array = np.array(coeff, dtype=object)
        result = np.dot(matrix, array)
        return(result.tolist())
    
    def pmap(self, coeff):
        return([x% abs(round((self.localization**2).real)) for x in coeff])
    
    def reduced(self, coeff, sde):
        if self.root_of_unity == 8:
            a = coeff[0] - coeff[4]
            b= coeff[1] - coeff[5]
            c = coeff[2] - coeff[6]
            d= coeff[3] - coeff[7]
            new_sde = sde
            while self.pmap([a,b,c,d]) == [0,0,0,0] and (a!= 0 or b!= 0 or c!= 0 or d!=0):
                a = round(a/2)
                b= round(b/2)
                c= round(c/2)
                d = round(d/2)
                new_sde += -2
            if self.pmap([a,b,c,d]) == [1,0,1,0] or self.pmap([a,b,c,d]) == [0,1,0,1] or self.pmap([a,b,c,d]) == [1,1,1,1]:
                return([round((b-d)/2), round((c+a)/2), round((b+d)/2), round((c-a)/2), 0,0,0,0], new_sde - 1)

            else:
                return([a,b,c,d,0,0,0,0], new_sde)

        else:
            reduced_coeff = coeff
            reduced_sde = sde
            while True and not(all(x == reduced_coeff[0] for x in reduced_coeff)):
                new_coeff = self.matrix(reduced_coeff, self.loc_char)
                if all(round(x)%(self.root_of_unity**2) == round(new_coeff[0])%(self.root_of_unity**2) for x in new_coeff):
                    reduced_coeff = [(round(x) - (round(new_coeff[0])%(self.root_of_unity**2)))//(self.root_of_unity**2) for x in new_coeff]
                    reduced_sde+=-1
                    
                else:
                    return(self.mode(reduced_coeff), reduced_sde)
            return(self.mode(coeff),sde)
        
    def mode(self, coeff):
        mode = max(set(coeff), key=lambda c: (coeff.count(c), c))
        return(self.add(coeff, [-mode] * self.num_coefficient))
    
    def subgroup(self, generators, depth = 10000):
        orbit = set()
        curr = random.choice(generators)
        for i in range(depth):
            curr = random.choice(generators) * curr
            orbit.add(curr)
        return(list(orbit))
    
    def subgroup_bfs(self, generators, depth = 10):
        orbit = set(generators)

        for i in range(depth):
            orbit = orbit.union(([g* o for g in generators for o in orbit]))
        
        return(list(orbit))

        

    def torus(self,subgroup, null_element):
        diags = []
        for elem in subgroup:
            if elem.is_diag(null_element):
                diags.append(elem)
        return(diags)
    
    def permutation_subgroup(self,subgroup, one_element, null_element=None):
        if null_element is None:
            null_element = cyclotomic_element(self, [0]*self.num_coefficient, 0)
        permuts = []
        for elem in subgroup:
            if elem.is_permutation(one_element, null_element):
                permuts.append(elem)
        return(permuts)

    
    def from_orbit(self, generator_set, depth= 100):
        curr = random.choice(generator_set)
        for i in range(depth):
            curr = random.choice(generator_set) * curr
        return(curr)
    
    def quotient(self,G, H, right = True):
        group = set(G)
        reps = []
        if right:
            while len(group)> 0 :
                new_elem = group.pop()
                coset = set([new_elem * h for h in H])

                group = group.difference(coset)
                reps.append(coset.pop())


            return(reps)
        else:
            while len(group)> 0 :
                new_elem = group.pop()
                coset = set([h*new_elem for h in H])

                group = group.difference(coset)
                reps.append(coset.pop())


            return(reps)

        




class cyclotomic_element:
    def __init__(self, ring, coefficients, sde = 0) -> None:
        self.ring = ring

        self.coefficients, self.sde = self.ring.reduced(coefficients,sde)

        if all([x==0 for x in self.coefficients]):
            self.sde = 0

    def __add__(self, value: object) -> object:

        if self.ring.root_of_unity == 8:
            denom = [0, 1, 0, 0, 0, 0, 0, 1]
        else:
            denom = gauss_sequence(self.ring.root_of_unity)

        if self.sde == value.sde:
            new_val = self.ring.add(self.coefficients, value.coefficients)
            return( cyclotomic_element(self.ring, new_val, self.sde))

        elif self.sde > value.sde:
            new_val = value.coefficients
            for i in range(self.sde - value.sde):
                new_val = self.ring.mul(new_val, denom)
            
            return(cyclotomic_element(self.ring, self.ring.add(new_val, self.coefficients), self.sde))
        else:
            new_self = self.coefficients
            for i in range(value.sde - self.sde):
                new_self = self.ring.mul(new_self, denom)
            
            return(cyclotomic_element(self.ring, self.ring.add(new_self, value.coefficients), value.sde))

    def __mul__(self, value: object) -> object:
        if type(value) == float or type(value) == int:
            if value == 0:
                return(cyclotomic_element(self.ring, [0] * self.ring.num_coefficient, self.sde))
            if type(value) == float and not value.is_integer():
                # Contract: non-integer scalars must be ±|λ|^k and are absorbed exactly into the sde.
                k = math.log(abs(value), abs(self.ring.localization))
                if not math.isclose(k, round(k)):
                    raise TypeError(f"cannot multiply ring element by {value}: not an integer or a power of the localization")
                k = round(k)
                sign = 1 if value > 0 else -1
                if self.ring.root_of_unity % 4 == 3:
                    if k % 2 != 0:
                        raise TypeError(f"cannot multiply ring element by {value}: odd powers of √p are not in the ring for p ≡ 3 (mod 4)")
                    if (k // 2) % 2 != 0:
                        sign = -sign
                return(cyclotomic_element(self.ring, [sign * c for c in self.coefficients], self.sde - k))
            return(cyclotomic_element(self.ring, [i * int(value) for i in self.coefficients], self.sde))
        else:
            return(cyclotomic_element(self.ring, self.ring.mul(self.coefficients, value.coefficients), self.sde + value.sde))
        
    def __rmul__(self, value):
        return(self*value)

    def power(self,value):
        result = self
        for i in range(value-1):
            result = result*self

        return(result)
    
    def conj(self):
        new_coeff = self.coefficients[1:]
        new_coeff.reverse()
        new_coeff = [self.coefficients[0]] + new_coeff
        # for p ≡ 3 (mod 4), conj(√p) = -√p, so odd sde picks up a sign flip
        if self.ring.root_of_unity % 4 == 3 and self.sde % 2 != 0:
            new_coeff = [-x for x in new_coeff]
        return(cyclotomic_element(self.ring, new_coeff, self.sde))
    
    def norm(self):
        return((self.comp() * (self.conj()).comp()).real)

    def comp(self):
        zeta = np.exp(2j * np.pi / self.ring.num_coefficient)
        return sum(c * zeta ** i for i, c in enumerate(self.coefficients)) / self.ring.localization ** self.sde

    def __eq__(self, other):
        if type(other) == cyclotomic_element:
            return(self.coefficients==other.coefficients and self.sde == other.sde)
        else:
            return(False)
        
    def is_monomial(self):
        nonzeros = []
        for x in self.coefficients:
            if x!=0:
                nonzeros.append(x)

        if len(nonzeros) == 1 or len(nonzeros) == 0:
            return(True)

        return(False)

    def __hash__(self):
        return hash((tuple(self.coefficients), self.sde))

    
    def __repr__(self):
        poly_string = ''
        for index in range(self.ring.num_coefficient):

            if self.coefficients[index] < 0:
                if index == 0:
                    poly_string += '-' + str(-self.coefficients[index])
                
                elif poly_string == '':
                    poly_string += '-' +  str(-self.coefficients[index]) + "\u03B6" + superscript_map.get(str(index))
                else:
                    poly_string += ' - ' +  str(-self.coefficients[index]) + "\u03B6" + superscript_map.get(str(index))
            elif self.coefficients[index] > 0:
                if index == 0:
                    poly_string +=  str(self.coefficients[index])
                elif poly_string == '':
                    poly_string += str(self.coefficients[index]) + "\u03B6" + superscript_map.get(str(index))
                else:
                    poly_string += ' + ' +  str(self.coefficients[index]) + "\u03B6" + superscript_map.get(str(index))


        return(poly_string)



class operator:
    def __init__(self, m, n, elements : Optional[list[list[cyclotomic_element]]]) -> None:
        if type(elements) == list:
            self.matrix = np.array([np.array(inner_list) for inner_list in elements])
        else:  
            self.matrix = elements
        self.m = m
        self.n = n
        self.sde = elements[0][0].sde
        self.shape = (m,n)

        self.string = ''

    def power(self, exponent):
        new_mat = self
        for i in range(exponent-1):
            new_mat = self*new_mat
        return(new_mat)




    def tensor(self, oper):
        return(operator(self.m * oper.m, self.n*oper.n, np.kron(self.matrix, oper.matrix)) )
    
    def tensor_power(self, power):
        final_matrix = self.matrix
        for i in range(power-1):
            final_matrix = np.kron(final_matrix, self.matrix)
        
        return(operator(self.m ** power, self.n**power, final_matrix))
            
        
    
    def __mul__(self, value):
        
        if type(value) == float or type(value) == int:
            return(operator(self.m, self.n, self.matrix * value))
        elif self.m == value.m and self.n ==1 and value.n ==1:
            ring = self.matrix[0][0].ring
            total = cyclotomic_element(ring, [0]*ring.num_coefficient, 0)
            for i in range(self.m):
                total = total + self.matrix[i][0].conj() * value.matrix[i][0]
            return(operator(1, 1, [[total]]))
        else:
            a = operator(self.m, value.n, np.matmul(self.matrix, value.matrix))
            a.string = self.string + value.string
            return(a)
        
    def __rmul__(self,value):
        if type(value) == float or type(value) == int:
            return(operator(self.m, self.n, self.matrix * value))
        
    def sde_profile(self):
        return(np.array([[obj.sde for obj in row] for row in self.matrix]))
    
    def sde_sum(self):
        return(np.sum(self.sde_profile()))

    
    def comp(self):
        return(np.array([[obj.comp() for obj in row] for row in self.matrix]))
    
    def unitary_check(self):
        res = np.dot(self.comp(), np.conjugate(self.comp().T))
        identity_matrix = np.eye(self.comp().shape[0])
        return(np.allclose(res, identity_matrix, atol=1e-8))
    
    def monomial_check(self):
        for x in self.matrix:
            for y in x:
                if not(y.is_monomial()):
                    return(False)
        if self.sde == 0:
            return(True)
        else:
            return(False)
    
    def synth_search(self,dropping_set):
        # first left-multiplier in dropping_set (ordered by priority) whose product has lower total sde
        for option in dropping_set:
            new_oper = option*self
            if new_oper < self:
                return(new_oper, option.string)

    def synthesize(self, dropping_set, target_sde = 1):
        mat = self
        final_string = ''
        while min(min(row) for row in mat.sde_profile()) > target_sde:
            result = mat.synth_search(dropping_set)
            if result is None:
                raise RuntimeError(f"synthesize: no dropping gate reduces SDE below {target_sde}")
            mat, string = result
            final_string = string + final_string

        return(final_string)





    def is_diag(self, null_element):
        for rows in range(self.m):
            for columns in range(self.n):
                if self.matrix[rows][columns] != null_element and rows != columns:
                    return(False)
        
        return(True)

    def is_permutation(self, one_element, null_element):

        for rows in range(self.m):
            one_counter = 0
            for columns in range(self.n):
                if self.matrix[rows][columns] == one_element:
                    one_counter+=1
                elif self.matrix[rows][columns] != null_element:
                    return(False)
            if one_counter != 1:
                return(False)
        
        return(True)



    def __lt__(self, other):
        return(self.sde_sum() < other.sde_sum())
    
    def __gt__(self, other):
        return(self.sde_sum() > other.sde_sum())
    
    def __eq__(self, other):
        return(np.array_equal(self.matrix, other.matrix))
    
    def __hash__(self):
        return hash(tuple(tuple(row) for row in self.matrix.tolist()))
        




    def __repr__(self):
        matrix = self

        rows = matrix.matrix.shape[0]
        placement  = rows//2 -1
        scalars = []
        for i in range(rows):
            if i == placement:
                scalars.append('√'+ str(round((matrix.matrix[0][0].ring.localization**2).real))+'^(-'+ str(matrix.sde) + ')')

            else:
                scalars.append('')
        headers = [''] + [f'Column {i}' for i in range(1, matrix.matrix.shape[1] + 1)]

        matrix_with_scalars = np.column_stack((scalars, matrix.matrix))

        return(tabulate(matrix_with_scalars, headers, tablefmt='fancy_grid'))

        
class state(operator):
    def __init__(self, d, unit_vector: list[cyclotomic_element]) -> None:
        rows = []
        for element in unit_vector:
            rows.append(np.array([element]))
        super().__init__(d, 1, np.array(rows))

    def norm(self):
        return((self*self).comp()[0][0])

