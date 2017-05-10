import numpy as np
import math as m
import matplotlib.pyplot as plt

imp0 = 377.0
size = 1800
# isOutput = False
isOutput = True
def fdtd_fresnel(epsilon1, epsilon2):
    n1 = np.sqrt(epsilon1)
    n2 = np.sqrt(epsilon2)
    eps = np.ones(size)
    eps[:int(size/2.0)] = epsilon1
    eps[int(size/2.0):] = epsilon2

    c = 1/np.sqrt(eps[0])
    al = (c-1)/(c+1)
    bl = 2/(c + 1)
    wl_nm1,wl_n,wl_np1 = 0,0,0 # Field at x=0 at time steps n-1, n, n+1
    wlp1_nm1,wlp1_n,wlp1_np1 = 0,0,0 # Field at x=1 at time steps n-1, n, n+1

    # Right boundary
    c = 1/np.sqrt(eps[-1])
    ar = (c-1)/(c+1)
    br = 2/(c + 1)
    wr_nm1,wr_n,wr_np1 = 0,0,0 # Field at x=size at time steps n-1, n, n+1
    wrm1_nm1,wrm1_n,wrm1_np1 = 0,0,0 # Field at x=size-1 at time steps n-1, n, n+1

    #Source
    source_width = 30.0*np.sqrt(max(epsilon1,epsilon2))
    delay = 10*source_width
    source_x = int(1.0*size/10.0)
    def source(current_time, delay, source_width):
        return m.exp(-(current_time-delay)**2/(2.0 * source_width**2))
    #Monitor points
    Emax1_x = int(2.0*size/15.0)
    Emax2_x = int(14.0*size/15.0)
    Emax1, Emax2 = 0,0
    Hmax1, Hmax2 = 0,0


    #Model
    total_steps = int(size*3.0+delay)  # Time stepping
    frame_interval = int(total_steps/15.0)
    all_steps = np.linspace(0, size-1, size)

    #Inital field E_z and H_y is equal to zero
    ez = np.zeros(size)
    hy = np.zeros(size)
    x = np.arange(0,size-1,1)
    #print(x)
    for t in range(total_steps):
        hy[x] = hy[x] + (ez[x+1] - ez[x])/imp0
        wrm1_np1 = hy[-2]
        wr_np1 = -wrm1_nm1 + ar*(wrm1_np1+wr_nm1) + br*(wr_n+wrm1_n)
        hy[-1] = wr_np1
        #Cycle field values at boundary
        wr_nm1, wrm1_nm1 = wr_n, wrm1_n
        wr_n, wrm1_n = wr_np1, wrm1_np1

        ez[x+1] = ez[x+1] + (hy[x+1]-hy[x])*imp0/eps[x+1]
        ez[source_x] += source(t, delay, source_width)

        #Evaluate Mur ABC value (eq. 6.35 Taflove)
        wlp1_np1 = ez[1]
        wl_np1 = -wlp1_nm1 + al*(wlp1_np1+wl_nm1) + bl*(wl_n+wlp1_n)
        ez[0] = wl_np1
        #Cycle field values at boundary
        wl_nm1, wlp1_nm1 = wl_n, wlp1_n
        wl_n, wlp1_n = wl_np1, wlp1_np1
        ######################
        #Monitor
        ######################
        Emax1 = max(Emax1, np.abs(ez[Emax1_x]))
        Emax2 = max(Emax2, np.abs(ez[Emax2_x]))
        Hmax1 = max(Hmax1, np.abs(hy[Emax1_x]))
        Hmax2 = max(Hmax2, np.abs(hy[Emax2_x]))
        ######################
        # Output
        ######################
        if t % frame_interval == 0 and isOutput:
            plt.plot(all_steps, ez, all_steps, hy*imp0)
            plt.pause(0.000001)
            plt.clf()

    print("n1 = %f, n2 = %f" % (n1, n2))
    Fresnel_ratio = 1-np.abs((n1-n2)/(n1+n2))**2
    print("Fresnel equation ratio 1-|(n1-n2)/(n1+n2)|^2 = %f" % Fresnel_ratio)
    FDTD_ratio = Emax2*Hmax2/(Emax1*Hmax1)
    print ("FDTD ratio = %f"%FDTD_ratio)
    error = np.abs((FDTD_ratio-Fresnel_ratio)/Fresnel_ratio)
    print("Error = %f%%" % (error*100.0) )
    return epsilon1, epsilon2, FDTD_ratio,Fresnel_ratio, error*100
epsilon1 = 8
epsilon2 = 2
a1=fdtd_fresnel(epsilon1, epsilon2)
epsilon1 = 3
epsilon2 = 5
a2=fdtd_fresnel(epsilon1, epsilon2)
epsilon1 = 1
epsilon2 = 9
a3=fdtd_fresnel(epsilon1, epsilon2)