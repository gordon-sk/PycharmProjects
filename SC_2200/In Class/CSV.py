import csv

x = open("lec26.csv", "rU")
a = csv.reader(x)


sheet = []
for row in a:
    sheet.append(row)
print()
print(sheet)

print()
sheet[3][3] = 100.00
print(sheet)
sum = 0
for x in range(1,4):
    sum += float(sheet[3][x])
sum /= 3
sheet[3][4] = round(sum,2)
print(sheet[3][4])