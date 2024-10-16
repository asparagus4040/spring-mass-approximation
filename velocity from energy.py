import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

# constant properties
m = 1
k = 1

# initial conditions
x0 = 0
v0 = 1

# simulation parameters
dt = 0.015
iterations = 100

# empty arrays
t = np.linspace(0, dt*(iterations-1), iterations)
x = np.zeros(iterations)
v = np.zeros(iterations)
a = np.zeros(iterations)
energy = np.zeros(iterations)

# array initial values
x[0] = x0
v[0] = v0
a[0] = -k/m * x0
energy[0] = 0.5*k*(x0**2) + 0.5*m*(v0**2)

for i in range(iterations - 1):
    # x(n+1) = x(n) + dt*v(n) + 0.5*dt**2*a(n)
    x[i+1] = x[i] + dt*v[i] + 0.5*(dt**2)*a[i]

    v[i+1] = sqrt( (2*energy[0]-k*x[i]**2)/m )

    # a(n+1) = -k/m * x(n+1)
    a[i+1] = -k/m * x[i+1]

    # energy should stay constant
    energy[i+1] = 0.5*k*(x[i+1]**2) + 0.5*m*(v[i+1]**2)

# plot data
fig, axs = plt.subplots(3, 2)

axs[0][0].plot(t, x, linestyle="solid", marker=".", color="orange")
axs[0][0].set_title("Position")
axs[0][0].set_xlabel("time (s)")
axs[0][0].set_ylabel("position (m)")

axs[0][1].plot(t, v, linestyle="solid", marker=".", color="red")
axs[0][1].set_title("Speed")
axs[0][1].set_xlabel("time (s)")
axs[0][1].set_ylabel("speed (m/s)")

axs[1][0].plot(t, a, linestyle="solid", marker=".", color="blue")
axs[1][0].set_title("Acceleration")
axs[1][0].set_xlabel("time (s)")
axs[1][0].set_ylabel("acceleration (m/s^2)")

axs[1][1].plot(t, energy, linestyle="solid", marker=".", color="green")
axs[1][1].set_title("Energy")
axs[1][1].set_xlabel("time (s)")
axs[1][1].set_ylabel("energy (J)")

axs[2][0].plot(x, v, linestyle="solid", marker=".", color="purple")
axs[2][0].set_title("Position and Speed")
axs[2][0].set_xlabel("position (m)")
axs[2][0].set_ylabel("speed (m/s)")

fig.suptitle("Spring-mass system")
fig.tight_layout()

plt.show()