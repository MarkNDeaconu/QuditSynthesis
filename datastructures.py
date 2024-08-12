from sde import *
import math 
import numpy as np
from typing import Optional

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
        self.loc_char = circulant(gauss_sequence(root_of_unity))

            
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
    
    def reduced(self, coeff):
        reduced_coeff = coeff
        while all(isinstance(x, int) for x in reduced_coeff):
            new_coeff = self.matrix(reduced_coeff, circulant(gauss_sequence(self.root_of_unity)))
            for i in range(len(reduced_coeff)):



class cyclotomic_element:
    def __init__(self, ring, coefficients, sde = 0) -> None:
        self.ring = ring
        self.coefficients, self.sde = reduced(coefficients,sde)

    def __add__(self, value: object) -> object:
        new_val = self.ring.add(self.coefficients, value.coefficients)

        return(cyclotomic_element(self.ring, new_val, self.sde))

    def __mul__(self, value: object) -> object:
        if type(value) == float or type(value) == int:
            scalar = value
            negative = False
            if value == 0:
                return(cyclotomic_element(self.ring, [0] * self.ring.num_coefficient, self.sde))
            elif value <0:
                negative = True
                scalar = -scalar
            
            if math.isclose(math.log(scalar, self.ring.localization), math.round(math.log(scalar, self.ring.localization))):
                new_sde = self.sde - math.round(math.log(value, self.ring.localization))
                if negative:
                    return(cyclotomic_element(self.ring, self.coefficients, new_sde))
                else:
                    return(cyclotomic_element(self.ring, self.ring.mul(self.coefficients, [-1] + [0]*(self.ring.num_coefficient-1)), new_sde))
            
            else:
                return cyclotomic_element(self.ring, [i * value for i in self.coefficients], self.sde)
        else:
            return(cyclotomic_element(self.ring, self.ring.mul(self.coefficients, value.coefficients), self.sde + value.sde))
        
    
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
    def __repr__(self):
        return(str(self.matrix))

        
class state(operator):
    def __init__(self, d, unit_vector: list[cyclotomic_element]) -> None:
        rows = []
        self.sde = unit_vector[0].sde
        for element in unit_vector:
            rows.append(np.array([element]))
        super().__init__(d, 1, np.array(rows))

    def norm(self):
        return(self*self)

