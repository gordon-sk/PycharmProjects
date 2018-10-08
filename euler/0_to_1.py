import random
import math

counts = []
for x in range(10000000):
    sum = 0
    count = 0
    while sum < 1:
        sum += random.random()
        count += 1
    counts.append(count)
    count = 0
    sum = 0

sum = 0
for x in counts:
    sum += x

print(sum / counts.__len__())
print(math.e)