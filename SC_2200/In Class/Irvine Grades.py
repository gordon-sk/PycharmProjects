import csv

def main():

    F0Bj = open("lec26.csv", "rU")
    read_data = csv.reader(F0Bj)
    for row in read_data:
        print(row)
    F0Bj.close()

    read_data[3][3] = '100.00'
    my_sum = 0





    pass


if __name__ == '__main__':
    main()