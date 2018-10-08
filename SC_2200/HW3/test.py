"""
 File name: test.py
 Author: Gordon Kiesling
 VUnetid: Kiesligs
 Email: Gordon.s.kiesling@vanderbilt.edu
 Class: CS2204
 Date: 2/26/16
 Honor statement: I have neither given nor recieved unauthorized aid on this work
 Assignment Number: 3
 Description: Unit tests of breast_cancer.py
     Use this file to test each of your functions individually.
"""

from breast_cancer import *


def file_io_test():
    the = get_diagnoses("results.txt")
    results = get_patient_data("patients.txt")
    pass


def gauss_test():
    # YOUR CODE HERE
    pass


def transpose_test():
    data = get_patient_data("patients.txt")
    results = get_diagnoses("results.txt")
    my_solve(data,results)
    pass

# The correctness of your my_solve function is part of the main analysis

if __name__ == "__main__":
    file_io_test()
    gauss_test()
    transpose_test()
