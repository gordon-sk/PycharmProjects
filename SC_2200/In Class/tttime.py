import time

def add(num):
    start = time.time()
    sum = 0
    for x in range(num + 1):
        sum += x
    end = time.time()

    return sum, round(end-start, 3)

if __name__ == "__main__":
    print(repr(add(7340000)) + " seconds" )