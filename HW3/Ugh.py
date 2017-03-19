import numpy
numpy.set_printoptions(linewidth=1000,suppress=True)

def transpose(matrix):
    """
    @param matrix - a numpy matrix
    @return transposed numpy matrix
    """
    new_Matrix = numpy.copy(matrix)
    return new_Matrix.transpose()

def gauss(a, b):
    """
    Note: be sure to read the project specification for this function
    @param a as matrix (2D numpy array)
    @param b as column vector (2D numpy array)
    @return result of the gauss as a numpy array
    """
    print("nothing has changed... let's gauss it")
    print(a)
    print(b)
    print()
    aug = numpy.zeros((a.shape[0],int(b.shape[0])+1))
    for x in range(a.shape[0]):
        for y in range(a.shape[0]):
            aug[x,y] = a[x,y]
    for z in range (b.shape[0]):
        aug[z,a.shape] = b[z]

    print()
    print("Nice, now we have our augmented matrix, lets solve")
    print(aug)
    print()
    print()
    # row 0 is the pivot row




    print("And now we have a solved matrix:")
    print(aug)
    return

def my_solve(a, b):
    """
    @param a training data as matrix (2D numpy array)
    @param b classifications as column vector (2D numpy array with only 1 column)
    @return solution as numpy array
    """
    A_T = transpose(a)
    A_T_dot_A = A_T.dot(a)
    A_T_dot_B = A_T.dot(b)

    print("a transposed")
    print(A_T)
    print()
    print("multiplied")
    print(A_T_dot_A)
    print(A_T_dot_B)
    print()

    solution_matrix = gauss(A_T_dot_A,A_T_dot_B)
    return solution_matrix

if __name__ == '__main__':

    a =  [[-3, 1], [1, 1], [-7, 1], [5 ,1]]
    b = [[70], [21], [110], [-35]]

    print("Initial values")
    print(a)
    print(b)
    print()

    w = my_solve(a,b)