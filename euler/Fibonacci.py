import os


def main():
    print(fib(10))

def fib(n):
    if n == 0:
        print("n is 0")
        return 0
    elif n == 1:
        print("n is 1")
        return 1
    else:
        print("n is " + str(n))
        return 0



main()