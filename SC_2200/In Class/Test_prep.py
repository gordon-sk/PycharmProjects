def palindrome_int(sample):
    new = str(sample)
    if is_palindrome(new):
        print("it's already a palindrome")
        print(sample)
    else:
        to_add = ''
        for x in range(1, new.__len__() + 1):
            to_add += new[-x]
        next = int(to_add) + int(new)
        is_palindrome(str(next))
        return next





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
    palindrome_int(88)
    palindrome_int(1342)