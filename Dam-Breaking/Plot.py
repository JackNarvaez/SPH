import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np

Boundary = np.loadtxt('./Files/Boundary.txt', dtype='float')
xb = Boundary[0:-2:2]
yb = Boundary[1:-1:2]

L = 4.0
H = 3.0

fig, ax = plt.subplots()
plt.scatter(xb, yb, color="r", zorder=3)
#plt.plot(xb, yb, marker='o', markersize=3, markeredgewidth=0., markerfacecolor='r')

Position = np.loadtxt('./Files/Position.txt', skiprows=0, max_rows=1, dtype='float')
Density = np.loadtxt('./Files/Density.txt', skiprows=0, max_rows=1, dtype='float')
x0 = Position[0:-2:2]
y0 = Position[1:-1:2]

sc = ax.scatter(x0, y0, color="b")

def plot(i):
    Position = np.loadtxt('./Files/Position.txt', skiprows=50*i, max_rows=1, dtype='float')
    Density = np.loadtxt('./Files/Density.txt', skiprows=50*i, max_rows=1, dtype='float')
    X = np.c_[Position[0:-2:2], Position[1:-1:2]]
    sc.set_offsets(X)
    #ax.plot(Position[0:-2:2], Position[1:-1:2], marker='o', color="b", ms=0.5)
    print(i)
    # manually relim:

anim = matplotlib.animation.FuncAnimation(fig, plot, frames=100, interval=10) 
anim.save('animation.mp4', fps=10)
#plt.show()
