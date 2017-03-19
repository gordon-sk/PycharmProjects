"""
# File name: Recursion.py
# Author: Gordon Kiesling
# VUnetid: Kiesligs
# Email: Gordon.s.kiesling@vanderbilt.edu
# Class: CS2204
# Assignment Number: 5
# Honor statement: I have neither given nor received unauthorized aid on this work
# Description: This set of functions perform a number of tasks using almost no lists.
#              Instead, recursion is used. Descriptions of each function are available
               within the code.
# Last Changed: 3/28/16

General note:
All NumPy array parameters and all list parameters will be
1-dimensional arrays/lists.
"""

import numpy


"""
Task: compute the sum of an array of numbers. You are required to use 
    a numpy array instead of the standard Python list for this exercise
Pre: arr is a 1-dimensional numpy array of integers
Post: the sum of arr[0]...arr[size-1] is returned
Additional requirement: You must do this by dividing the array in half on 
each recursive call rather than reducing the array size by only one
"""


def sum_array(arr):
    if arr.size == 0:
        return 0
    elif arr.size == 1:
        return arr[0]
    else:
        return sum_array(arr[:arr.size/2]) + sum_array((arr[arr.size/2:arr.size]))

"""
Task: determine if target is in the set
Pre: my_set is a list representing a set of objects
Post: true is returned if target is in my_set, else false;
      mySet is unchanged
"""


def member(target, my_list):
    if my_list.__len__() == 0:
        foo = False
    else:
        if my_list[my_list.__len__()-1] == target:
            foo = True
        else:
            foo = member(target, my_list[:my_list.__len__()-1])
    return foo


"""
Task: compute the sum of the first n harmonic terms
Pre: n is a positive integer
Post: the sum of the first n harmonic terms is returned.
    The harmonic series is 1 + (1/2) + (1/3) + (1/4) + ...
"""


def harmonic_sum(n):
    count = 0
    if n == 0:
        pass
    else:
        if str(n).__contains__("9"):
            count += 0 + harmonic_sum(n)
        else:
            count += 1/n + harmonic_sum(n-1)
    return count

"""
Task: determine if a string is a palindrome
Pre: my_str is a string object
Post: returns true if my_str is a palindrome, otherwise returns false
      The test is case insensitive [look up upper() and lower() of the 
      string class]. You do not need to worry about
      trimming blanks from the ends of the string.
"""


def is_palindrome(my_str):
    if my_str.__len__() == 0:
        foo = True
    else:
        if my_str[0].upper() == my_str[-1].upper():
            foo = is_palindrome(my_str[1:-1])
        else:
            foo = False
    return foo


"""
Task: replace all occurrences of 'target' in the array 'numbers'
      with 'replacement'
Pre: 'numbers' is a NumPy array of integers
Post: all occurrences of 'target' in 'numbers' have been replaced
      with 'replacement'; the number of replacements performed is
      returned to the caller.
"""


def replace(target, replacement, numbers, count=0):
    if numbers.size == 0:
        pass
    else:
        if numbers[-1] == target:
            count += 1
            numbers[-1] = replacement
        count = replace(target, replacement, numbers[:-1], count)
    return count


"""
Task: compute the Greatest Common Divisor (GCD) of two non-negative
        integers using Euclid's formula:
 
Euclid's method for computing the greatest common divisor (GCD) of two 
non-negative integers a and b is as follows. Divide a and b to obtain the
integer quotient q and remainder r, so that a = bq+r (if b = 0, 
then GCD(a, b) = a). Then GCD(a, b) = GCD(b, r). Replace a with b and 
b with r and repeat the procedure. Because the remainders are decreasing, 
eventually a remainder of 0 will result. The last nonzero remainder is 
the greatest common divisor of a and b.
 
Pre: the parameters x & y are non-negative
Post: the GCD of x & y is returned
"""


def gcd(x, y):
    to_return = 0
    b = min(x, y)
    a = max(x, y)
    if a < 0 or b < 0:
        print("Error: Negative value(s) passed")
    elif b == 0:
        to_return = a
    else:
        q = int(a/b)
        r = a % b
        if q == 0:
            to_return = a
        else:
            to_return = gcd(b, r)
    return to_return

"""
 Task: reverse the contents of an array of numbers 
     You are to use a NumPy array instead of the standard Python list or this exercise
 Pre: 'arr' is a (potentially empty) NumPy array of integers
 Post: the elements of the array have been reversed.
"""


def reverse_list(arr):
    if arr.__len__() == 0:
        pass
    else:
        plchold = arr[0]
        arr[0] = arr[-1]
        arr[-1] = plchold
        reverse_list(arr[1:-1])


"""
 Task: produce the binary representation of a decimal number
   A decimal number is converted to binary by repeated
   division by 2. For each division, keep track of the quotient
   and remainder. The remainder becomes the low-order bit (rightmost
   bit) of the binary representation. The higher-order bits are
   determined by repeating the processes with the quotient.
   The process stops when num is either zero or one.
 Pre: num is a non-negative integer
 Post: the binary representation of num is produced and returned
   as a string.    
"""


def convert2binary(num):
    binary = ""
    if num == 0 or num == 1:
        binary += repr(num) + binary
    else:
        x = int(num/2)
        r = int(num % 2)
        binary += str(convert2binary(x)) + repr(r)
    return binary

    
"""
 Task: Print a pseudo hourglass pattern on the screen
 Pre: num is a positive integer
 Post: the desired pattern is displayed on the screen
 Use Python string operations to print a single line of *'s of
 a desired size. Use recursion to complete the pattern.
 Example: a call to print_pattern(4) should produce the
 7-line pattern below:

****
***
**
*
**
***
****

"""


def print_pattern(num, plc=0):

    if plc == 0:
        plc = num
    if num == 1:
        print_pattern(num-2, plc)
    elif num > 1:
        for x in range(num):
            print("*", end="")
        print()
        print_pattern(num-1, plc)
        pass
    elif (-1*plc) <= num <= -1:
        num *= -1
        for x in range(num):
            print("*", end="")
        print()
        print_pattern((-1*num)-1, plc)

"""
 Task: initialize all elements of the array between indices lb and ub to the
   given value, including the elements at lb & ub
 Note: lb = lower bound, ub = upper bound
     You are to use the NumPy array instead of the standard Python list for this exercise
 Pre: lb and ub are valid indices into the NumPy array arr [the actual size of the array is
   not important]
 Post: the array elements in the segment arr[lb..ub] have been initialized to 'value'
 Additional requirement: This function must be done by dividing the array segment
   in half and performing recursive calls on each half (as opposed to just shrinking
   the array bound by one each time)
"""


def array_initialize(arr, value, lb, ub):
    arr = arr[lb:ub]
    if arr.size == 0:
        pass
    else:
        arr[0] = value
        array_initialize(arr[1:], value, 0, arr.size)


"""
 Task: Compute the Binomial Coefficient using Pascal's Triangle.
 The Binomial Coefficient B(n, r) is the coefficient of the term x^r in the binomial
 expansion of (1 + x)^n. For example, B(4, 2) = 6 because (1+x)^4 = 1 + 4x + 6x^2 + 4x^3 + x^4. 
 The Binomial Coefficient can be computed using the Pascal Triangle formula:
    B(n, r) = 1                          if r = 0 or r = n
    B(n, r) = B(n-1, r-1) + B(n-1, r)    otherwise
 Pre: r & n are non-negative, and r<=n
"""


def binomial_coeff(n, r):
    if n < 0 or r < 0:
        print("error, r or n is negative")
    elif r == 0 or r == n:
        return 1
        pass
    else:
        return binomial_coeff(n-1, r-1) + binomial_coeff(n-1, r)

"""
 Task: Given a string, compute recursively a new string where all the adjacent 
  chars are now separated by a "*".
 Pre: my_str is a string (may be empty).
 Post: a correctly starred string is returned.
 Examples:  
   add_star("hello") --> "h*e*l*l*o"
   add_star("abc") --> "a*b*c"
   add_star("ab") --> "a*b"
"""


def add_star(my_str):
    new_str = ''
    if my_str.__len__() == 1:
        new_str += my_str[-1]
    elif my_str.__len__() != 0:
        new_str += my_str[0] + "*" + add_star(my_str[1:])
    return new_str


"""
 Task: Given a non-negative int n, compute recursively (no loops) the count of the 
   occurrences of 2 as a digit, except that a 2 with another 2 immediately to its 
   left counts double, so 2212 yields 4. Note that mod (%) by 10 yields the rightmost 
   digit (126 % 10 is 6), while divide (/) by 10 removes the rightmost digit (126 / 10 is 12). 
 Pre: n is non-negative
 Post: count of the occurrences of 2 is returned (with doubling as appropriate)
 Examples:
   count2(2) --> 1
   count2(212) --> 2
   count2(2212) --> 4
"""


def count2(n, count=0):
    if n == 0:
        pass
    else:
        if n % 10 == 2:
            count += 1
            if (int(n/10)) % 10 == 2:
                count += 1
        count = count2(int(n/10), count)
    return count


"""
 Task: Given a string and a non-empty substring sub, compute recursively the number 
   of times that sub appears in the string, without the sub strings overlapping. 
 Pre: sub is a non-empty string
 Post: the count of non-overlapping occurrences of sub in str is returned
 Examples:
   count_subs("catcowcat", "cat") --> 2
   count_subs("catcowcat", "cow") --> 1
   count_subs("catcowcat", "dog") --> 0
"""


def count_subs(my_str, sub, count=0):
    if my_str.__len__() == 0:
        pass
    elif my_str[0] == sub[0]:
        if my_str[:sub.__len__()] == sub:
            count += 1
        count = count_subs(my_str[sub.__len__():], sub, count)
    else:
        count = count_subs(my_str[1:], sub, count)
    return count


"""
 Task: Given a string, compute recursively a new string where all the lowercase 'x' chars 
   have been moved to the end of the string. 
 Pre: my_str is a string (possibly empty)
 Post: return a new string where all lowercase 'x' chars have been moved to the end 
 Examples:
   move_xs("xxre") --> "rexx"
   move_xs("xxhixx") --> "hixxxx"
   move_xs("xhixhix") --> "hihixxx"
"""


def move_xs(my_str, x_str='', new_str=''):
    if my_str.__len__() == 0:
        pass
    else:
        if my_str[0] == 'x':
            x_str += 'x'
        else:
            new_str += my_str[0]
        new_str = move_xs(my_str[1:], '', new_str) + x_str
    return new_str


"""
Our main method for testing all the recursive functions above.
It is your job to add tests so that your code is fully tested.
"""


def main():

    print("Adding contents of list: ")
    arr1 = numpy.array([1, 2, 3, 4, 5])
    print("Summation of [1,2,3,4,5] = ", end="")
    print(sum_array(arr1))
    arr2 = numpy.array([6, 7, 8, 9, 10, 11, 12])
    print("Summation of [6,7,8,9,10,11,12] = ", end="")
    print(sum_array(arr2))
    arr3 = numpy.array([6, -7, 8, 9, -10, 11, 12])
    print("Summation of [6,-7,8,9,-10,11,12] = ", end="")
    print(sum_array(arr3))
    arr4 = numpy.array([])
    print("Summation of an empty array = ", end="")
    print(sum_array(arr4))
    print()

    print("Finding a target in a list of numbers: ")
    print("4 in set [1,2,3,4,5]?: ", end="")
    print(member(4, [1, 2, 3, 4, 5]))
    print("6 in set [1,2,3,4,5]?: ", end="")
    print(member(6, [1, 2, 3, 4, 5]))
    print("10 in an empty set?: ", end="")
    print(member(10, []))
    print()

    print("Calculating harmonic sum: ")
    print("Sum of 5 is: ", end="")
    print(harmonic_sum(5))
    print("Sum of 100 is: ", end="")
    print(harmonic_sum(100))
    print("Sum of 950 is: ", end="")
    print(harmonic_sum(950))
    print()

    print("Checking for palindrome: ")
    print("America? ", end="")
    print(is_palindrome("America"))
    print("madam? ", end="")
    print(is_palindrome("madam"))
    print("deleveled? ", end="")
    print(is_palindrome("deleveled"))
    print("Hey there amigo? ", end="")
    print(is_palindrome("Hey there amigo"))
    print("Empty string? ", end="")
    print(is_palindrome(''))
    print()

    print("Replacing numbers from array: ")
    arr = numpy.array([2, 3, 4, 5, 6, 7, 5])
    print("Replacing 5 with 7 in " + repr(arr) + " produces ", end="")
    print(repr(replace(5, 7, arr)) + " replacements. \nArray is now: ", end="")
    print(arr)
    print()
    arr1 = numpy.array([2, 3, 4, 5, 6, 7, 5])
    print("Replacing 10 with 7 in " + repr(arr1) + " produces ", end="")
    print(repr(replace(10, 7, arr1)) + " replacements. \nArray is now: ", end="")
    print(arr1)
    print()

    print("Finding GCD: ")
    print("GCD of 24, 36 is: ", end="")  # should be 12
    print(gcd(24, 36))
    print("GCD of 112, 378 is: ", end="")  # should be 14
    print(gcd(112, 378))
    print()

    print("Reverse list:")
    lst = numpy.array([1, 2, 3, 4, 5, 6, 7])
    print(repr(lst) + " reversed is: ", end="")
    reverse_list(lst)
    print(lst)
    print()

    print("Convert to binary:")
    print("Binary rep of 14 is: ", end="")
    print(convert2binary(14))
    print("Binary rep of 1234 is: ", end="")
    print(convert2binary(1234))
    print()

    print("Check printPattern: ")
    print_pattern(4)
    print()
    print()

    print("Testing arrayInitialize: ")
    arr = numpy.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    print("Indices 2 through 6 initialized to 0 in " + repr(arr2) + "\nlooks like ", end='')
    array_initialize(arr2, 0, 2, 6)
    print(arr)
    print()

    print("Checking binomial coefficient:")
    print("Binomial coefficient of 4,2 is: ", end="")
    print(binomial_coeff(4, 2))
    print("Binomial coefficient of 5,3 is: ", end="")
    print(binomial_coeff(5, 3))
    print("Binomial coefficient of 50,3 is: ", end="")
    print(binomial_coeff(50, 3), end="")
    print(" (That's actually correct!)")
    print()

    print("Checking addStar: ")
    print("'welcome' converts to: ", end="")
    print(add_star("welcome"))
    print("'ABCDE' converts to: ", end="")
    print(add_star("ABCDE"))
    print()

    print("Checking count2: ")
    print("2212 through count2 equals: ", end="")
    print(count2(2212))
    print("2221212 through count2 equals: ", end="")
    print(count2(2221212))
    print()

    print("Checking countSubs: ")
    print("Number of times 'is'' appears in 'This is a fun fishing dish': ", end="")
    print(count_subs("This is a fun fishing dish", "is"))
    print("Number of times 'cow'' appears in 'catcowcowcatcow': ", end="")
    print(count_subs("catcowcowcatcow", "cow"))
    print()

    print("Checking moveXs: ")
    print("xhixhixx becomes: ", end="")
    print(move_xs("xhixhixx"))
    print("'Ix've gxotx xa bxad feexlixng' becomes: ", end="")
    print(move_xs("Ix've gxotx xa bxad feexlixng"))
    print()



    print()
    print()
    print()
    print("Check this shit.")
    print(harmonic_sum(100))


""" If we are the executable, call the main() method """
if __name__ == '__main__':
    main()