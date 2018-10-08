import mpl
num = 600851475143
sqrt = m.sqrt(num)
potentials = []
for x in range(int(sqrt+1)):
    if x == 0 or x == 1 or x % 2 == 0:
        pass
    else:
        if num % x == 0:
            potentials.append(x)

print(potentials)
prime = True
finals = []
for x in potentials:
    for y in range(2, x):
        if x%y == 0:
            prime = False
            break
    if prime:
      finals.append(x)
print(finals)