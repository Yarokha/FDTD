import numpy as np
import math as m
import matplotlib.pyplot as plt

imp0 = 377.0
size = 1800
epsilon = 5

c = 1/np.sqrt(epsilon)
a = (c-1)/(c+1)
b = 2/(c + 1)
# Left boundary
wl_nm1,wl_n,wl_np1 = 0,0,0 # Field at x=0 at time steps n-1, n, n+1
wlp1_nm1,wlp1_n,wlp1_np1 = 0,0,0 # Field at x=1 at time steps n-1, n, n+1
# Right boundary
wr_nm1,wr_n,wr_np1 = 0,0,0 # Field at x=size at time steps n-1, n, n+1
wrm1_nm1,wrm1_n,wrm1_np1 = 0,0,0 # Field at x=size-1 at time steps n-1, n, n+1

source_width = 45.0*np.sqrt(epsilon)
delay = 10*source_width
source_x = int(1.0*size/2.0)

def source(current_time, delay, source_width):
    return m.exp(-(current_time-delay)**2/(2.0 * source_width**2))

total_steps = int(1*(size+delay)*np.sqrt(epsilon))
frame_interval = int(total_steps/15.0)
all_steps = np.linspace(0, size-1, size)
ez = np.zeros(size)
hy = np.zeros(size)
x = np.arange(0, size-1, 1)

for t in range(total_steps):
    hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0
    wrm1_np1 = hy[-2]
    wr_np1 = -wrm1_nm1 + a*(wrm1_np1+wr_nm1) + b*(wr_n+wrm1_n)
    hy[-1] = wr_np1
    wr_nm1, wrm1_nm1 = wr_n, wrm1_n
    wr_n, wrm1_n = wr_np1, wrm1_np1
    ez[x+1] = ez[x+1] + (hy[x+1]-hy[x])*imp0/epsilon
    ez[source_x] += source(t, delay, source_width)
    #Evaluate Mur ABC value (eq. 6.35 Taflove)
    wlp1_np1 = ez[1]
    wl_np1 = -wlp1_nm1 + a*(wlp1_np1+wl_nm1) + b*(wl_n+wlp1_n)
    ez[0] = wl_np1
    #Cycle field values at boundary
    wl_nm1, wlp1_nm1 = wl_n, wlp1_n
    wl_n, wlp1_n = wl_np1, wlp1_np1
    if t % 100 == 0:
        plt.plot(all_steps , ez, all_steps, hy*imp0)
        plt.pause(0.00001)
        plt.clf()
