number = open("13_number.txt")
data = number.read()
number.close()

data = data.splitlines()
sum = 0
for x in data:
    sum += int(x)
print(str(sum)[:10])

print("Now we are running a different script from the command line")