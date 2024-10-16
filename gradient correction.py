import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

# constant properties
m = 1
k = 3

# initial conditions
x0 = 1
v0 = 1

# simulation parameters
dt = 0.1
iterations = 200
correction_factor = 0.01
max_correction_steps = 20

# empty arrays
t = np.linspace(0, dt*(iterations-1), iterations)
x = np.zeros(iterations)
v = np.zeros(iterations)
a = np.zeros(iterations)
energy = np.zeros(iterations)
correction_steps = np.zeros(iterations)

# array initial values
x[0] = x0
v[0] = v0
a[0] = -k/m * x0
energy[0] = 0.5*k*(x0**2) + 0.5*m*(v0**2)

for i in range(iterations - 1):
    x[i+1] = x[i] + dt*v[i] + 0.5*(dt**2)*a[i]

    v[i+1] = v[i] + dt*a[i]

    # energy should stay constant
    energy[i+1] = 0.5*k*(x[i+1]**2) + 0.5*m*(v[i+1]**2)

    # correction (only goes downwards)
    current_cor_steps = 0
    while energy[i+1] > energy[0] and current_cor_steps < max_correction_steps:
        gradient_magnitude = sqrt( (k*x[i+1])**2 + (m*v[i+1])**2 )

        x[i+1] = x[i+1] - correction_factor * k * x[i+1] / gradient_magnitude
        v[i+1] = v[i+1] - correction_factor * m * v[i+1] / gradient_magnitude

        energy[i+1] = 0.5*k*(x[i+1]**2) + 0.5*m*(v[i+1]**2)
        current_cor_steps += 1
    
    # store number of correction steps needed
    correction_steps[i+1] = current_cor_steps

    a[i+1] = -k/m * x[i+1]

# get max energy deviation
energy_deviation = (max(energy) - min(energy))/energy[0]
print(energy_deviation*100, "%")

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

axs[2][1].plot(t, correction_steps, linestyle="solid", marker=".", color="magenta")
axs[2][1].set_title("Number of correction steps required")
axs[2][1].set_xlabel("time (s)")
axs[2][1].set_ylabel("correction steps")

fig.suptitle("Spring-mass system")
fig.tight_layout()

plt.show()