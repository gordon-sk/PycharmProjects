# Author: Gordon Kiesling

import lec15_dllist as llist

if __name__ == '__main__':

    A = llist.DLList()

    for x in range(8):
        A.append(x)
    print(A)

    for x in range(8):
        str = input("Please provide a key and value: ")
