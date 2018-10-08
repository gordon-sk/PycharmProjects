f1 = 0
f2 = 1
fib_list = [f1, f2]

while str(fib_list[-1]).__len__() < 3:
    fib_list.append(fib_list[-1] + fib_list[-2])

for x in fib_list:
    print()
    print(x)
print(fib_list.__len__())