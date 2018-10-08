# Gordon Kiesling
# SC3250
# HW6 -- Magnetic Pendulum

import matplotlib.pyplot as plt
import scipy.integrate as intg
import numpy as np
import math

magnets = []
magnets.append([-.250,  .433])
magnets.append([-.250, -.433])
magnets.append([.500,  .00])
R, C, d = .20, .50, .25
x0, y0 = .2, .2
vx0, vy0 = 0, 0
init_values = [vx0, vy0, x0, y0] 

# Defining timegrid for odeint
maxT = 100
timestep = .01
nTimeSteps = int(maxT/timestep)
timeGrid = np.linspace(0, maxT, nTimeSteps)

# Defining the derivative function
def derivative(var_list, t):
    vx = var_list[0]
    vy = var_list[1]
    x = var_list[2]
    y = var_list[3]
    
    dvxdt = -1 * R * vx - C * x
    dvydt = -1 * R * vy - C * y
    for k in magnets:
        top1 = k[0] - x
        top2 = k[1] - y
        bottom = ((k[0]-x) ** 2 + (k[1]-y) ** 2 + d ** 2) ** 1.5
        dvxdt += top1/bottom
        dvydt += top2/bottom
    return [dvxdt, dvydt, vx, vy]



# Finding Scipy's solution (for reference)
odeint_sol = intg.odeint(derivative, [0, 0, -.2, -.2], timeGrid)
x, y = [], []
for l in odeint_sol:
    x.append(l[2])
    y.append(l[3])
plt.plot(x, y, color='orange')
plt.plot(magnets[0][0], magnets[0][1], 'or', zorder=10, label='magnet 1')
plt.plot(magnets[1][0], magnets[1][1], 'ob', zorder=10, label='magnet 2')
plt.plot(magnets[2][0], magnets[2][1], 'og', zorder=10, label='magnet 3')
plt.legend()
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.grid(True)
plt.xlim(-.3, .55)
plt.ylim(-.5, .5)
plt.title("Scipy odeint solution")
plt.show()



# In[3]:

# Part 1 -- RK4 Solution
# Prof told me to 'do it all out' because my list comprehension was whack. I'm not apologizing 
# for this obscene codeblock. Sorry :(

print("hi")

def magnet_runge_kutta(f, inits, t):
    value_grid = [[0, 0, 0, 0] for k in range(t.__len__())]
    value_grid[0] = inits # [vx, vy, x, y] for n iterations
    for n in range(len(value_grid) - 1):
        # x value
        X1 = timestep *  f([value_grid[n][0], value_grid[n][1], value_grid[n][2], value_grid[n][3]], t[n])[2]
        Y1 = timestep *  f([value_grid[n][0], value_grid[n][1], value_grid[n][2], value_grid[n][3]], t[n])[3]
        VX1 = timestep * f([value_grid[n][0], value_grid[n][1], value_grid[n][2], value_grid[n][3]], t[n])[0]
        VY1 = timestep * f([value_grid[n][0], value_grid[n][1], value_grid[n][2], value_grid[n][3]], t[n])[1]
        
        X2 =  timestep * f([value_grid[n][0] + .5 * VX1, value_grid[n][1] + .5 * VY1, 
                            value_grid[n][2] + .5 * X1,  value_grid[n][3] + .5 * Y1], t[n] + .5 * timestep)[2]
        Y2 =  timestep * f([value_grid[n][0] + .5 * VX1, value_grid[n][1] + .5 * VY1, 
                            value_grid[n][2] + .5 * X1,  value_grid[n][3] + .5 * Y1], t[n]+ .5 * timestep)[3]
        VX2 = timestep * f([value_grid[n][0] + .5 * VX1, value_grid[n][1] + .5 * VY1, 
                            value_grid[n][2] + .5 * X1,  value_grid[n][3] + .5 * Y1], t[n]+ .5 * timestep)[0]
        VY2 = timestep * f([value_grid[n][0] + .5 * VX1, value_grid[n][1] + .5 * VY1, 
                            value_grid[n][2] + .5 * X1,  value_grid[n][3] + .5 * Y1], t[n]+ .5 * timestep)[1]
        
        X3 =  timestep * f([value_grid[n][0] + .5 * VX2, value_grid[n][1] + .5 * VY2, 
                            value_grid[n][2] + .5 * X2,  value_grid[n][3] + .5 * Y2], t[n]+ .5 * timestep)[2]
        Y3 =  timestep * f([value_grid[n][0] + .5 * VX2, value_grid[n][1] + .5 * VY2, 
                            value_grid[n][2] + .5 * X2,  value_grid[n][3] + .5 * Y2], t[n]+ .5 * timestep)[3]
        VX3 = timestep * f([value_grid[n][0] + .5 * VX2, value_grid[n][1] + .5 * VY2, 
                            value_grid[n][2] + .5 * X2,  value_grid[n][3] + .5 * Y2], t[n]+ .5 * timestep)[0]
        VY3 = timestep * f([value_grid[n][0] + .5 * VX2, value_grid[n][1] + .5 * VY2, 
                            value_grid[n][2] + .5 * X2,  value_grid[n][3] + .5 * Y2], t[n]+ .5 * timestep)[1]
        
        X4 =  timestep * f([value_grid[n][0] + VX3, value_grid[n][1] + VY3, 
                            value_grid[n][2] + X3,  value_grid[n][3] + Y3], t[n] + timestep)[2]
        Y4 =  timestep * f([value_grid[n][0] + VX3, value_grid[n][1] + VY3, 
                            value_grid[n][2] + X3,  value_grid[n][3] + Y3], t[n] + timestep)[3]
        VX4 = timestep * f([value_grid[n][0] + VX3, value_grid[n][1] + VY3, 
                            value_grid[n][2] + X3,  value_grid[n][3] + Y3], t[n] + timestep)[0]
        VY4 = timestep * f([value_grid[n][0] + VX3, value_grid[n][1] + VY3, 
                            value_grid[n][2] + X3,  value_grid[n][3] + Y3], t[n] + timestep)[1]
        
        x  = value_grid[n][2] + ((1/6) * (X1 + 2 * X2 + 2 * X3 + X4))
        y  = value_grid[n][3] + ((1/6) * (Y1 + 2 * Y2 + 2 * Y3 + Y4))
        vx = value_grid[n][0] + ((1/6) * (VX1 + 2 * VX2 + 2 * VX3 + VX4))
        vy = value_grid[n][1] + ((1/6) * (VY1 + 2 * VY2 + 2 * VY3 + VY4))
        value_grid[n+1] = [vx, vy, x, y]
        
    return value_grid


# In[4]:

# Part 2
# Returns INT index of settled-upon magnet. Actual magnet label in plots is +1 of this number
def settle_magnet(x, y):
    radii = []
    for k in magnets:
        r = math.sqrt((k[0] - x) ** 2 + (k[1] - y)**2)
        radii.append(r)
    return radii.index(min(radii))


# In[7]:

# Part 3 -- color coded plotting

color_dict = {0:'r', 1:'b', 2:'g'}
x1, y1, x2, y2, x3, y3 = [], [], [], [], [], []

init_values1 = [0, 0, .2, .2]
rk4_sol1 = magnet_runge_kutta(derivative, init_values1, timeGrid)
for l in rk4_sol1:
    x1.append(l[2])
    y1.append(l[3])

init_values2 = [0, 0, -.2, -.2]
rk4_sol2 = magnet_runge_kutta(derivative, init_values2, timeGrid)
for l in rk4_sol2:
    x2.append(l[2])
    y2.append(l[3])
    
init_values3 = [0, 0, 0, .4]
rk4_sol3 = magnet_runge_kutta(derivative, init_values3, timeGrid)
for l in rk4_sol3:
    x3.append(l[2])
    y3.append(l[3])

    
plt.plot(x1, y1, color=color_dict[settle_magnet(rk4_sol1[-1][2], rk4_sol1[-1][3])], label='solution 1')
plt.plot(x2, y2, color=color_dict[settle_magnet(rk4_sol2[-1][2], rk4_sol2[-1][3])], label='solution 2')
plt.plot(x3, y3, color=color_dict[settle_magnet(rk4_sol3[-1][2], rk4_sol3[-1][3])], label='solution 3')
plt.annotate(xy=(magnets[0][0], magnets[0][1]), zorder=10, label='magnet 1', 
             ha='center', va='center', s='+', color='red')
plt.annotate(xy=(magnets[1][0], magnets[1][1]), zorder=10, label='magnet 2', 
             ha='center', va='center', s='+', color='blue')
plt.annotate(xy=(magnets[2][0], magnets[2][1]), zorder=10, label='magnet 3', 
             ha='center', va='center', s='+', color='green')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.legend()
plt.grid(True)
plt.title("Custom RK4 Solution")
plt.show()

print("You can't actually see the magnets because the assignment called for color coding\n",
      "them the same color as the paths.\n")
print('starting positions are:')
print('[vx0, vy0, x0, y0 : solution x]')
print(init_values1, ': solution 1, magnet 3')
print(init_values2, ': solution 2, magnet 2')
print(init_values3, ': solution 3, magnet 1')


# In[15]:

# Part 4 pt 1

init_values = [0, 0, -.2, -.2]
rk4_4_1 = magnet_runge_kutta(derivative, init_values, timeGrid)
k1, i1 = [], []
for x in rk4_new:
    k1.append(x[2])
    i1.append(x[3])
'''
init_values = [0, 0, 0, 0]
rk4_4_2 = magnet_runge_kutta(derivative, init_values, timeGrid)
k2, i2 = [], []
for x in rk4_new:
    k2.append(x[2])
    i2.append(x[3])
init_values = [0, 0, .2, .2]
rk4_4_3 = magnet_runge_kutta(derivative, init_values, timeGrid)
k3, i3 = [], []
for x in rk4_new:
    k3.append(x[2])
    i3.append(x[3])
'''
plt.plot(k1, i1)
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.legend()
plt.grid(True)
plt.annotate(xy=(magnets[0][0], magnets[0][1]), zorder=10, label='magnet 1', 
             ha='center', va='center', s='+', color='red')
plt.annotate(xy=(magnets[1][0], magnets[1][1]), zorder=10, label='magnet 2', 
             ha='center', va='center', s='+', color='blue')
plt.annotate(xy=(magnets[2][0], magnets[2][1]), zorder=10, label='magnet 3', 
             ha='center', va='center', s='+', color='green')
plt.title("Custom RK4 Solution")
plt.show()


# In[ ]:



