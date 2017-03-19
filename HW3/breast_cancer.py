"""
 File name: breast_cancer.py
 Author: Gordon Kiesling
 VUnetid: Kiesligs
 Email: Gordon.s.kiesling@vanderbilt.edu
 Class: CS2204
 Date: 1/26/16
 Honor statement: I have neither given nor received unauthorized aid on this work
 Assignment Number: 3
 Description: Given a set of training data, generates a
     prediction vector for breast cancer biopsies.
"""

import numpy
import timeit
numpy.set_printoptions(linewidth=1000, suppress=True)

"""
Be sure to read thoroughly the project specification as it
provides necessary information pertinent to the understanding
and implementation of these methods.

While Python allows for object oriented design, it is also a
very popular scripting language. This program will differ from 
your previous two in that it will not make use of the class structure.
"""


def get_patient_data(filename):
    """
    @param filename as string
    @return numpy array of data
    """
    fin = open(filename)
    count = 0
    for row in fin:
        count += 1          # counting number of rows
    total_list = [[] for x in range(count)]
    fin.close()

    data_string = ""
    fin = open(filename)
    for row in fin:
        data_string += row
    fin.close()

    list_data = data_string.splitlines()    # necessary to avoid /n char
    for x in range(list_data.__len__()):
        total_list[x] = list_data[x].split(" ")

    for y in range(total_list.__len__()-1):
        for z in range(10):
            total_list[y][z] = float(total_list[y][z])
    data_array = numpy.zeros((count, 10))

    for a in range(total_list.__len__()):
        for b in range(10):
                data_array[a, b] = total_list[a][b]

    return data_array

def get_diagnoses(filename):
    """
    @param filename as string
    @return numpy array of data
    """
    data_string = ""
    fin = open(filename)
    for file_row in fin:
        data_string += file_row
    fin.close()
    return numpy.array(data_string.splitlines(), dtype='float64')


def transpose(matrix):
    """
    @param matrix - a numpy matrix
    @return transposed numpy matrix
    """
    new_matrix = numpy.copy(matrix)
    return new_matrix.transpose()

def gauss(a, b):
    """
    Note: be sure to read the project specification for this function
    @param a as matrix (2D numpy array)
    @param b as column vector (2D numpy array)
    @return result of the gauss as a numpy array
    """
    aug = numpy.zeros((a.shape[0],int(b.shape[0])+1))
    for x in range(a.shape[0]):
        for y in range(a.shape[0]):
            aug[x, y] = a[x, y]
    for z in range (b.shape[0]):
        aug[z, a.shape] = b[z]
    print(aug)

    # row 0 is the pivot row
    for pivot_row in range(10):
        for row_ in range(10):
            pivot_value = aug[pivot_row, pivot_row]
            if row_ == pivot_row:
                for x in range(11):
                    aug[pivot_row, x] /= pivot_value
            else:
                for column_ in range(pivot_row, 11):
                    scaling_value = aug[row_, pivot_row] / pivot_value
                    aug[row_, column_] -= scaling_value*aug[pivot_row, column_]
    print(aug)
    solution_array = numpy.zeros((10))
    for add in range (10):
        solution_array[add] = aug[add,10]
    return solution_array

def my_solve(a, b):
    """
    @param a training data as matrix (2D numpy array)
    @param b classifications as column vector (2D numpy array with only 1 column)
    @return solution as numpy array
    """
    A_T = transpose(a)
    A_T_dot_A = A_T.dot(a)
    A_T_dot_B = A_T.dot(b)

    solution_matrix = gauss(A_T_dot_A, A_T_dot_B)
    return solution_matrix


# normally we would have this line call our main method, but to make
# calls to timeit simpler, we will implement our driver code below
if __name__ == '__main__':
    #a =  [[-3, 1], [1, 1], [-7, 1], [5 ,1]]
    #b = [[70], [21], [110], [-35]]
    # solution to the above equation should be
    # c = [[-12.1],[29.4]]

    # load the training set
    a = get_patient_data('patients.txt')
    b = get_diagnoses('results.txt')

    # first, time the execution of your my_solve vs numpy's solution
    x = timeit.timeit('numpy.linalg.lstsq(a,b)', 'from __main__ import a,b; import numpy', number=1)
    y = timeit.timeit('my_solve(a,b)', 'from __main__ import a,b,my_solve', number=1)
    print('numpy:', x)
    print('custom:', y)

    # solve both ways
    # z is the learned knowledge base / decision set
    z = my_solve(a, b)
    # w is the solution from numpy.linalg
    w = numpy.linalg.lstsq(a, b)[0]
    print()
    print("Results compared...")
    print(z)
    print("Z")
    print(w)
    print("W")

    # test your solution against the numpy solution
    # Percent correct (if < 1.0 (100%) you need to check your work)
    correct = 0
    if z.shape != w.shape:
        print('Regression performed incorrectly.')
    for i in range(z.shape[0]):
        if abs(z[i]-w[i]) < .00001:
            correct += 1

    print('\nLibrary Test')
    print('percent correct: ', float(correct) / 10 )

    # test your code against the test set
    # load the test set
    a = get_patient_data("patients1.txt")
    b = get_diagnoses("results1.txt")

    # Check the percent correct
    # Percent incorrect
    # Percent of total that are false positives
    # Percent of total that are false negatives
    correct = 0
    f_pos = 0
    f_neg = 0
    for row in range(a.shape[0]):
        sol = 0
        for col in range(a.shape[1]):
            sol += a[row][col]*z[col][0]
        sol = 0 if sol < 0.5 else 1
        if sol == b[row][0]:
            correct += 1
        elif sol == 1:
            f_pos += 1
        else:
            f_neg += 1

    print('\nData Test')
    print('percent correct: ', (float(correct)/a.shape[0]))
    print('false positives: ', (float(f_pos)/a.shape[0]))
    print('false negatives: ', (float(f_neg)/a.shape[0]))
