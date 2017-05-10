import numpy as np
import math as m
import matplotlib.pyplot as plt

imp0 = 377.0
size = 1800

epsilon = 5

source_width = 30.0*np.sqrt(epsilon)
delay = 10*source_width
source_x = int(1.0*size/2.0)

def source(current_time, delay, source_width):
    return m.exp(-(current_time-delay)**2/(2.0 * source_width**2))

total_steps = int((size+delay)*np.sqrt(epsilon))
frame_interval = int(total_steps/30.0)
all_steps = np.linspace(0, size-1, size)

ez = np.zeros(size)
hy = np.zeros(size)
x = np.arange(0, size-1, 1)
for t in range(total_steps):
    hy[-1] = hy[-2]
    hy[x] = hy[x] + (ez[x+1] - ez[x]) / imp0
    ez[0] = ez[1]
    ez[x + 1] = ez[x + 1] + (hy[x + 1] - hy[x]) * imp0 / epsilon
    ez[source_x] += source(t, delay, source_width)
    if t % 100 == 0:
        plt.plot(all_steps, ez)
        plt.plot(all_steps, hy*imp0)
        plt.pause(0.00001)
        plt.clf()
