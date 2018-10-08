# Also by Stew
import math


def is_palindrome(my_str):
    if my_str.__len__() == 0:
        foo = True
    else:
        if my_str[0].upper() == my_str[-1].upper():
            foo = is_palindrome(my_str[1:-1])
        else:
            foo = False
    return foo


def check_prime(n):
    if n <= 1:
        return False
    elif n <= 3:
        return True
    elif n%2 == 0 or n%3 == 0:
        return False
    i = 5
    while i*i <= n:
        if n%i == 0 or n%(i+2) == 0:
            return False
        i += 6
    return True

def count_divisors(n):
    count = 0
    for x in range(1, int(math.sqrt(n)) + 1):
        if n%x == 0:
            count += 2
    return count


def str_to_int(my_str):
    if my_str.isdigit():
        result = 0
        for x in range(my_str.__len__()):
            result += (float(my_str[x]) - 48) * math.pow(10, my_str.__len__()-x-1)
    return int(result)




