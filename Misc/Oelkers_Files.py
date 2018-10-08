import os
import matplotlib.pyplot as plt
filepath = '/Users/Stew/Desktop/Fall17/Light_curves/'
file_list = os.listdir(filepath)

for file in file_list:
    if file[-4:] == '.dat':
        f = open(filepath + file)
        data = f.read()
        f.close()
        data = data.splitlines()
        x, y = [], []
        for line in data:
            k = line.split(' ')
            if k.__contains__(''):
                k.remove('')
            date = float(k[0].strip())
            brightness = float(k[1].strip())
            x.append(date)
            y.append(brightness)
        plt.scatter(x, y, s=5)
        plt.grid(True)
        plt.gca().invert_yaxis()
        plt.title(" ")
        plt.savefig(filepath + file[:-4] + "_graph.png")
        plt.clf()