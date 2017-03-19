def three():
    num = 13195
    for x in range(num):
        prime = True
        for y in range(1,x):
            if x%y == 0:

                break
def double():
    x = 1
    while x <= 600851475143:
        x *= 2
        k = str(x)



def is_palindrome(my_str):
    if my_str.__len__() == 0:
        foo = True
    else:
        if my_str[0].upper() == my_str[-1].upper():
            foo = is_palindrome(my_str[1:-1])
        else:
            foo = False
    return foo


if __name__ == "__main__":
    double()