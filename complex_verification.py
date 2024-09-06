import numpy as np

def verify(string, dictionary, expected_answer=0):
    for key, value in dictionary.items():
        locals()[key] = value
    return(eval(string))