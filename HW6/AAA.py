fname_dict = {}
file = open("texts/hamlet.txt")
data = file.read()
file.close()
file_list = data.splitlines()

print(file_list)
print()
print(file_list[0].split())