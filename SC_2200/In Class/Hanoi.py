def Hanoi(n, a, b, c):
    if n == 1:
        print("Move top disk from", a, "to", c)
        return 1
    else:
        count = Hanoi(n-1, a, b, c)
        count += Hanoi(1, a, c, b)
        count += Hanoi(n-1, b, c, a)
        return count


if __name__ == '__main__':

    n = 10

    a = [10,9,8,7,6,5,4,3,2,1]
    print(a)
    b = []
    c = []
    Hanoi(n,a,b,c,)