import numpy as np
from matplotlib import pyplot as plt
import time
plt.ion()

def Ax(V,mask):
    #Computes A*x based on a potential and a mask
    Vuse=V.copy()
    Vuse[mask]=0
    ans=(Vuse[1:-1,:-2]+Vuse[1:-1,2:]+Vuse[2:,1:-1]+Vuse[:-2,1:-1])/4.0
    ans=ans-V[1:-1,1:-1]
    return ans

def pad(A):
    AA=np.zeros([A.shape[0]+2,A.shape[1]+2])
    AA[1:-1,1:-1]=A
    return AA

def cg(n, tol):
    V=np.zeros([n,n])
    bc=0*V

    #Create mask for walls
    mask=np.zeros([n,n],dtype='bool')
    mask[:,0]=True
    mask[:,-1]=True
    mask[0,:]=True
    mask[-1,:]=True

    #Modify mask to also include the indices of a circle of radius r
    r=n//8
    x=np.arange(0,n)
    y=np.arange(0,n)
    cx=len(x)//2
    cy=len(y)//2
    circle_indices=(x[np.newaxis,:]-cx)**2+(y[:,np.newaxis]-cy)**2 < r**2
    bc[circle_indices]=1
    mask[circle_indices]=True

    #Compute b
    b=-(bc[1:-1,0:-2]+bc[1:-1,2:]+bc[:-2,1:-1]+bc[2:,1:-1])/4.0
    V=0*bc
    r=b-Ax(V,mask)
    p=r.copy()
    converged=False
    steps=0

    #While the residuals are smaller thn the tolerance, apply the conjugate gradient method to
    #try to solve for potential
    while(not converged):
        Ap=(Ax(pad(p),mask))
        rtr=np.sum(r*r)
        alpha=rtr/np.sum(Ap*p)
        V=V+pad(alpha*p)
        rnew=r-alpha*Ap
        beta=np.sum(rnew*rnew)/rtr
        p=rnew+beta*p
        r=rnew
        steps += 1
        converged = np.sum(r*r)<tol
        print('On iteration ', steps, ', res are ', np.sum(r*r))
    return V, steps
    
#Set a size of box and a tolerance
n=512
tol=0.01

#Keep track of time and apply the conjugate gradient method
start_time = time.time()
V, steps = cg(n,tol)
rho=V[1:-1,1:-1]-(V[1:-1,0:-2]+V[1:-1,2:]+V[:-2,1:-1]+V[2:,1:-1])/4.0
print('The conjugate gradient took ' + str(time.time() - start_time) + 's to run, corresponding to ' + str(steps) + ' steps.')

#Plot the final V and charge density
plt.figure()
plt.imshow(V)
plt.colorbar()
plt.title('Numerical V with conjugate gradient')
plt.savefig('q2_final_V.pdf')
plt.figure()
plt.imshow(rho)
plt.colorbar()
plt.title('Final charge distribution')
plt.savefig('q2_final_rho.pdf')