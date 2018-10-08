from scipy import special as S
import matplotlib.pyplot as plt

# part 1
x_vals, y_vals = [], []
x = 0
while x < 25:
    y_vals.append(S.y0(x))
    x_vals.append(x)
    x += .1

root = 0
for i in range(1, y_vals.__len__()):
    if y_vals[i-1] * y_vals[i] > 0:
        root = (x_vals[i], y_vals[i])
print("ugh")
plt.plot(x_vals, y_vals)
plt.grid(True)
plt.xlabel("x-values")
plt.ylabel("y-values")
plt.savefig("ugh.png")