def main():
    Julians = []
    Julians.append(2)
    ya = prompt()
    print("The year is", ya)

def prompt():
    year = input("What year: ")
    try:
        print("trying")
        year = int(year)
        if year < 1600 or year > 2200:
            print("out of range\n")
            year = prompt()
            return year
        else:
            print('nice')
            return year
    except ValueError:
        print("excepted")
        print("Not an integer\n")
        year = prompt()
        return year

main()