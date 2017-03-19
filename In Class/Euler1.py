# By Stew
import math
import Helper_methods_euler


def one():
    our_sum = 0
    for a in range(1000):
        if a % 3 == 0 or a % 5 == 0:
            our_sum += a
    return our_sum


def two():
    x = 0
    y = 1
    _sum = 0
    while _sum < 4000000:
        if y % 2 == 0:
            _sum += y
        plc = y
        y += x
        x = plc
    return _sum


def three():
    n = 600851475143
    x = n
    while x > 0:
        if n % x == 0:
            if Helper_methods_euler.check_prime(n):
                break
        x -= 1
    print(n)


def four():
    x = 999
    y = 999
    l = []

    for a in range(1, x):
        for b in range(1, y):
            if Helper_methods_euler.is_palindrome(str(a*b)):
                print(str(a*b) + "    " + str(a) + "    " + str(b))
                l.append(a*b)
        print()
        print()
        print()
    l.sort()
    print(l)


def five():
    n = 20
    x = 11
    found = False

    while not found:
        a = True
        for attempt in range(1, n):
            if x % attempt != 0:
                a = False
        if a:
            print(x)
            found = True
        else:
            x += 1


def six():
    n = 100
    sum1 = 0
    sum2 = 0

    for x in range(1, n+1):
        sum1 += math.pow(x, 2)
    for y in range(1, n+1):
        sum2 += y
    sum2 = math.pow(sum2, 2)
    print(sum2-sum1)


def seven():
    l = []
    n = 2
    while l.__len__() < 10005:
        if Helper_methods_euler.check_prime(n):
            l.append(n)
        n += 1
    print(l[10000])


def eight():
    file = open("8_data.txt")
    data = file.read()
    file.close()
    l = data.splitlines()
    num = ""
    for a in range(l.__len__()):
        num += l[a]

    results = []
    n = 13
    for b in range(num.__len__()-n):
        temp = num[b:b+n]
        print(temp)
        prod = 1
        for c in range(temp.__len__()):
            prod *= int(temp[c])
        results.append(prod)
    results.sort()
    print(results)


def nine():
    l = []
    for x in range(1, 500):
        print(repr(x) + " of 1000")
        for y in range(1, 500):
            for z in range(1, 500):
                if math.pow(x, 2) + math.pow(y, 2) == math.pow(z, 2) and x + y + z == 1000:
                    l.append([x, y, z])
    print(l)


def ten():
    _sum = 0
    for x in range(int(2*math.pow(10, 6))):
        if Helper_methods_euler.check_prime(x):
            _sum += x
        if x % 100 == 0:
            print(repr(x) + " down...")
    print("Done!")
    print(_sum)


def twelve():
    n = 1
    found = False
    while not found:
        n += 1
        delta = 0
        for x in range(n):
            delta += x
        y = Helper_methods_euler.count_divisors(delta)
        print(y)
        if y > 500:
            found = True
    print(delta)


def thirteen():
    file = open("prob_13.txt")
    bulk = file.read()
    file.close()
    data = bulk.splitlines()
    for x in range(data.__len__()):
        data[x] = float(data[x])
    print(str(data[0]))

    _sum = 0
    for x in range(10):
        _sum += data[x]
    print(_sum)


def fourteen():
    record = []
    max_count = 0
    for n in range(2, 1000000):
        count = 0
        plc = n
        while n != 1:
            if n % 2 == 0:
                n /= 2
            else:
                n = 3*n + 1
            count += 1
        if count > max_count:
            record = [count, plc]
            max_count = count
    print(record)


def fifteen():
    pass


def sixteen():
    base = math.pow(2, 1000)
    base = str(base)
    print(base)
    _sum = 0
    for x in range(base.__len__()-2):
        _sum += int(base[x])
    print(_sum)


def seventeen():
    n = 5
    storage = []
    for x in range(n):
        x += 1






if __name__ == '__main__':
    sixteen()