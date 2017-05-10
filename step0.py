import numpy as np
import matplotlib.pyplot as plt
import math as m

def checker(e,n):
    for i in range(n):
        if abs(e[i]) < 10**(-4):
            ch = False
        else:
            return True
    return ch

N = 500
sigma = 10.0
MaxTime = 1000
delay = 100.0
hy = np.zeros(N)
ez = np.zeros(N)
imp0 = 377.0
x = np.linspace(0, N-1, N)
# plt.ion()


for t in range(MaxTime):
    hy[:-1] += (ez[1:] - ez[:-1]) / imp0
    ez[1:] += (hy[1:] - hy[:-1]) * imp0
    ez[int(N/2)-1] += m.exp(-((t-delay)/sigma)**2)

    if t % 25 == 0:
        plt.plot(x, ez)
        plt.plot(x, hy*imp0)
        plt.pause(0.00001)
        plt.clf()
    # if not checker(ez, N) and t > 80000:
    #     plt.clf()
    #     plt.subplot(211)
    #     plt.plot(x, ez)
    #     plt.subplot(212)
    #     plt.plot(x, hy * imp0)
    #     plt.show()
    #     break
