list1 = [1,2,3,4,5,6,7,8,9]

def binary_search(set, target):
    search_location = set.__len__()
    x = 0
    while set[search_location] != target:
        x = x + 2
        if set[search_location] > target:
            search_location = search_location / 2
        elif set[search_location] < target:
            search_location = search_location + (x*set.__len__())


    return search_location

target = 3
answer = binary_search(list1,target)
print(answer)