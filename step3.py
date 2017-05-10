import numpy as np
import matplotlib.pyplot as plt

imp0 = 377.0
size = 1800
epsilon = 5

source_width = 30.0*np.sqrt(epsilon)
delay = 10*source_width

source_x = int(size/2.0)
def source(current_time, delay, source_width):
    return np.exp(-(current_time-delay)**2/(2.0 * source_width**2))

total_steps = int((size*2.0+delay)*np.sqrt(epsilon))
frame_interval = int(total_steps/15.0)
all_steps = np.linspace(0, size-1, size)

ez = np.zeros(size)
hy = np.zeros(size)

dx = 1
R_0 = 1e-5
m = 3.85
pml_width = 50
sxmax = -(m+1)*np.log(R_0)/2/imp0/(pml_width*dx)
sx = np.zeros(size)
sxm = np.zeros(size)
psi_hy = np.zeros(size)
psi_ez = np.zeros(size)

for i in range(int(pml_width)):
    sx[i+1] = sxmax*((pml_width-i-0.5)/pml_width)**m
    sxm[i] = sxmax*((pml_width-i)/pml_width)**m
    sx[size-i-1] = sxmax*((pml_width-i-0.5)/pml_width)**m
    sxm[size-i-1] = sxmax*((pml_width-i)/pml_width)**m

aez = np.exp(-sx*imp0)-1
bez = np.exp(-sx*imp0)
ahy = np.exp(-sxm*imp0)-1
bhy = np.exp(-sxm*imp0)

x = np.arange(0, size-1, 1)

for time in range(total_steps+1):
    psi_hy[x] = bhy[x]*psi_hy[x] + ahy[x]*(ez[x+1] - ez[x])
    hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0 + psi_hy[x]/imp0
    psi_ez[x+1] = bez[x+1]*psi_ez[x+1] + aez[x+1]*(hy[x+1]-hy[x])
    ez[x+1] = ez[x+1] + (hy[x+1]-hy[x])*imp0/epsilon + psi_ez[x+1]*imp0/epsilon
    ez[source_x] += source(time, delay, source_width)

    if time % frame_interval == 0:
        plt.plot(all_steps, ez, all_steps, hy*imp0)
        plt.pause(0.000001)
        plt.clf()
