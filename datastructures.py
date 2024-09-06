import math 
import numpy as np
from typing import Optional
from tabulate import tabulate

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


class cyclotomic_ring:
    def __init__(self, root_of_unity, localization) -> None:
        self.root_of_unity = root_of_unity
        self.localization = localization
        self.num_coefficient = root_of_unity
        self.loc_char = root_of_unity**2*np.linalg.inv(circulant(gauss_sequence(root_of_unity)))

            
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
        array = np.array(coeff)
        result = np.dot(matrix, array)
        return(result.tolist())
    
    def pmap(self, coeff):
        return([x% round((self.localization**2).real) for x in coeff[0:4]])
    
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
            else:
                return(self.mode(coeff),sde)
        
    def mode(self, coeff):
        mode = max(set(coeff), key=coeff.count)
        return(self.add(coeff, [-mode] * self.num_coefficient))



class cyclotomic_element:
    def __init__(self, ring, coefficients, sde = 0) -> None:
        self.ring = ring

        self.coefficients, self.sde = self.ring.reduced(coefficients,sde)

        if all([x==0 for x in coefficients]):
            self.sde = 0
        # self.coefficients, self.sde = coefficients,sde

    def __add__(self, value: object) -> object:
        if self.sde == value.sde:
            new_val = self.ring.add(self.coefficients, value.coefficients)
            return( cyclotomic_element(self.ring, new_val, self.sde))
        elif self.sde > value.sde:
            new_val = value.coefficients
            for i in range(self.sde - value.sde):
                new_val = self.ring.mul(new_val, gauss_sequence(self.ring.root_of_unity))
            
            return(cyclotomic_element(self.ring, self.ring.add(new_val, self.coefficients), self.sde))
        else:
            new_self = self.coefficients
            for i in range(value.sde - self.sde):
                new_self = self.ring.mul(new_self, gauss_sequence(self.ring.root_of_unity))
            
            return(cyclotomic_element(self.ring, self.ring.add(new_self, value.coefficients), value.sde))

    def __mul__(self, value: object) -> object:
        if type(value) == float or type(value) == int:
            scalar = value
            negative = False
            if value == 0:
                return(cyclotomic_element(self.ring, [0] * self.ring.num_coefficient, self.sde))
            elif value <0:
                negative = True
                scalar = -scalar

            
            
            if math.isclose(math.log(scalar, abs(self.ring.localization)), round(math.log(scalar, abs(self.ring.localization)))):
                new_sde = self.sde - round(math.log(scalar, abs(self.ring.localization)))
                if not(negative):
                    return(cyclotomic_element(self.ring, self.coefficients, new_sde))
                else:
                    return(cyclotomic_element(self.ring, self.ring.mul(self.coefficients, [-1] + [0]*(self.ring.num_coefficient-1)), new_sde))
            
            else:
                return cyclotomic_element(self.ring, [i * value for i in self.coefficients], self.sde)
        else:
            return(cyclotomic_element(self.ring, self.ring.mul(self.coefficients, value.coefficients), self.sde + value.sde))
        
    def __rmul__(self, value):
        return(self*value)

    def power(self,value):
        result = self
        for i in range(value-1):
            result = result*self

        return(result)

    def comp(self):
        zeta = np.exp(2j * np.pi / self.ring.num_coefficient)

        complex_code = np.diag([zeta**n for n in range(self.ring.num_coefficient)])

        return(sum(self.ring.matrix(self.coefficients, complex_code)))
    
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
        try:
            self.sde2 = elements[0][1].sde
        except:
            self.sde2 = 0
        self.shape = (m,n)

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
            left_vector = self.matrix[:, 0]
            right_vector = value.matrix[:, 0]
            return(float(np.vdot(left_vector, right_vector).complex()))
        else:
            return(operator(self.m, value.n, np.matmul(self.matrix, value.matrix)))
        
    def __rmul__(self,value):
        if type(value) == float or type(value) == int:
            return(operator(self.m, self.n, self.matrix * value))
        
    def sde_profile(self):
        return(np.array([[obj.sde for obj in row] for row in self.matrix]))
    
    def comp(self):
        return(np.array([[obj.comp()/(obj.ring.localization**obj.sde) for obj in row] for row in self.matrix]))

    def __repr__(self):
        matrix = self
        if type(matrix) == int or type(matrix) == float:
            return(matrix)
        else:

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
        self.sde = unit_vector[0].sde
        for element in unit_vector:
            rows.append(np.array([element]))
        super().__init__(d, 1, np.array(rows))

    def norm(self):
        return(self*self)


