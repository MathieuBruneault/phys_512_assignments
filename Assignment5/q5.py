import numpy as np
import matplotlib.pyplot as plt

#Set some parameters and variables that will be useful
dt = 0.0005
dy = 0.0005
k = 10**(-4)
y_max = 0.04
t_max = 5   
T0 = 100

def FTCS(dt,dy,t_max,y_max,k,T0):
    #Solve numerically the 1-D heat equation using Newmann B.C on the right wall
    #and linearly increasing temperature on the left wall
    s = k*dt/dy**2
    y = np.arange(0,y_max+dy,dy) 
    t = np.arange(0,t_max+dt,dt)
    r = len(t)
    c = len(y)
    T = np.zeros([r,c])
    #This is to increase linearly T on the left wall
    T[:,0] = np.arange(0,r)*dt
    #Loop throuh time and space
    for n in range(0,r-1):
        for j in range(1,c-1):
            T[n+1,j] = T[n,j] + s*(T[n,j-1] - 2*T[n,j] + T[n,j+1]) 
        #Newmann B.C for the right wall
        j = c-1 
        T[n+1, j] = T[n,j] + s*(T[n,j-1] - 2*T[n,j] + T[n,j-1])
    return y,T,r,s

#Get results
y,T,r,s = FTCS(dt,dy,t_max,y_max,k,T0)

#Loop through each time step ad plot
plot_times = np.arange(0.01,t_max,0.01)
for t in plot_times:
    plt.clf()
    plt.plot(y,T[int(t//dt),:])
    plt.ylim([0,np.amax(T[:,0])])
    plt.title('Temperature distribution')
    plt.xlabel('x')
    plt.ylabel('T')
    if (t==1 or t==2 or t==3):
        plt.savefig('q5_temp_distribution_at_t_' + str(t) + '.pdf')
    plt.pause(0.001)