def factorial(n):
    sum = 1
    for x in range(1, n+1):
        sum *= x
    return sum

k = factorial(100)
k = str(k)

finalsum = 0
for x in k:
    finalsum += int(x)

print(finalsum)
