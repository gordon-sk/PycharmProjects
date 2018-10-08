num = 2 ** 1000
num = str(num)
sum = 0
for x in range(num.__len__()):
    sum += int(num[x])
print(sum)