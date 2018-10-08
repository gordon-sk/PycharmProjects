import string
import random

def monkeys():
    a = "methinks it is like a weasel"
    shake = string.ascii_letters + " "
    new = ""
    timer = 0
    for y in range(1000):
        alist = {}
        while not a == new:
            timer += 1
            new = ''
            for x in range(27):
                new += random.choice(shake)
            score = 0
            for x in range(new.__len__()):
                if  new[x] == a[x]:
                    score += 1
            alist[new] = score
        print(alist)

if __name__ == "__main__":
    monkeys()