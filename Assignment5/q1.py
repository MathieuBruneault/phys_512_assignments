import numpy as np
from matplotlib import pyplot as plt
import time
plt.ion()

def Ax_2d(V,mask):
    #Compute A*x for some potential and mask
    Vuse=V.copy()
    Vuse[mask]=0
    ans=(Vuse[1:-1,:-2]+Vuse[1:-1,2:]+Vuse[2:,1:-1]+Vuse[:-2,1:-1])/4.0
    ans=ans-V[1:-1,1:-1]
    return ans

def relaxation_method(n,tol):
    #Get initial values for potential and set boundary conditions on the edges
    V=np.zeros([n,n])
    bc=0*V
    mask=np.zeros([n,n],dtype='bool')
    mask[:,0]=True
    mask[:,-1]=True
    mask[0,:]=True
    mask[-1,:]=True

    #Get the indices of a circle in the center of our space, set its potential to be fixed at 1
    r=n//8
    x=np.arange(0,n)
    y=np.arange(0,n)
    cx=len(x)//2
    cy=len(y)//2
    circle_indices=(x[np.newaxis,:]-cx)**2+(y[:,np.newaxis]-cy)**2 < r**2
    bc[circle_indices]=1
    mask[circle_indices]=True

    #Set some initial values
    V=bc.copy()
    b=-(bc[1:-1,0:-2]+bc[1:-1,2:]+bc[:-2,1:-1]+bc[2:,1:-1])/4.0
    converged=False
    steps=0

    while(not converged):
        #While the residuals are too big, compute the potential and new residuals using the elaxation method
        V[1:-1,1:-1]=(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0
        V[mask]=bc[mask]
        steps+=1
        r=b-Ax_2d(V,mask)
        converged=np.sum(r*r)<tol
        if steps%1000 == 0:
            print('On iteration ', steps, ', res are ', np.sum(r*r))
    return V, steps


#Set some initial values
n=512
tol=0.01
start_time = time.time()
V, steps = relaxation_method(n,tol)
rho=V[1:-1,1:-1]-(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0
print('The relaxation solver took ' + str(time.time() - start_time) + 's to run, corresponding to ' + str(steps) + ' steps.')

#Plot the final potential as well as the charge density
plt.figure()
plt.imshow(V)
plt.colorbar()
plt.savefig('q1_final_V.pdf')
plt.figure()
plt.imshow(rho)
plt.colorbar()
plt.savefig('q1_final_rho.pdf')