import random
set = [1,4,2]

for x in range(0,100):
    set = set + [random.randint(0,100)]

def bubble_sort(data):
    placeholder = 0

    for y in range(0,set.__len__()):
        for x in range (0,set.__len__()-1):
            if set[x] > set[x+1]:
                placeholder = set[x+1]
                set[x+1] = set[x]
                set[x] = placeholder

print(set)
bubble_sort(set)
print("Sorting...")
print(set)