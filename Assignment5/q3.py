import numpy as np
from matplotlib import pyplot as plt
import time

def make_rhs_2d(pot,mask):
    """Make the right hand side.  We know that in charge-free regions,
    V0 equals the average of its neighbors, or V0+0.25*(V_l+V_r +V_u+V_d)=0.  If
    some of the neighbors have been set by the boundary conditions, we have
    V0 + 0.25*sum(V interior) = -0.25*sum(V_boundary) and we only solve for 
    the potential not on the specified boundary conditions.  Note that if we had charge in there, 
    then up to a multiplicative constant, V0+0.25*(V_l+V_r +V_u+V_d)=rho, so you could add the 
    (suitably scaled) charge to the rhs."""
    mat=np.zeros(pot.shape)
    mat[:,:-1]=mat[:,:-1]+pot[:,1:]
    mat[:,1:]=mat[:,1:]+pot[:,:-1]
    mat[:-1,:]=mat[:-1,:]+pot[1:,:]
    mat[1:,:]=mat[1:,:]+pot[:-1,:]
    mat[mask]=0  #for the masked region, i.e. specified by the boundary conditions, we'll peg the RHS to zero since we keep 
                 #the potential fixed anyways
    return mat

def ax_2d(mat,mask,copy=False):
    """Write the Laplacian operator in the way we need it to be.  Note that the boundary conditions as specified by the mask
    do not enter into the matrix since they are on the right-hand side of the matrix equation.  So, set them to zero here, then we won't have
    to worry about handling them separately."""
    if copy:
        mat=mat.copy()
    mat[mask]=0
    mm=4*mat
    mm[:,:-1]=mm[:,:-1]-mat[:,1:]
    mm[:,1:]=mm[:,1:]-mat[:,:-1]
    mm[1:,:]=mm[1:,:]-mat[:-1,:]
    mm[:-1,:]=mm[:-1,:]-mat[1:,:]
    mm[mask]=0
    return mm
    
def cg(rhs,x0,mask,tol):
    """cg(rhs,x0,mask,niter) - this runs a conjugate gradient solver to solve Ax=b where A
    is the Laplacian operator interpreted as a matrix, and b is the contribution from the 
    boundary conditions.  Incidentally, one could add charge into the region by adding it
    to b (the right-hand side or rhs variable)"""

    t1=time.time()
    Ax=ax_2d(x0,mask)
    r=rhs-Ax
    p=r.copy()
    x=x0.copy()
    rsqr=np.sum(r*r)
    print('starting rsqr is ',rsqr)
    converged = False
    steps=0
    while(not converged):
        Ap=ax_2d(p,mask)
        alpha=np.sum(r*r)/np.sum(Ap*p)
        x=x+alpha*p
        r=r-alpha*Ap
        rsqr_new=np.sum(r*r)
        beta=rsqr_new/rsqr
        p=r+beta*p
        rsqr=rsqr_new
        steps+=1
        converged=rsqr<tol
    t2=time.time()
    print('final rsqr is ',rsqr,' after ',t2-t1,' seconds')
    return x

def deres_mat(mat):
    """A quick and dirty way to downgrade the resolution of a potential matrix by a 
    factor of 2.  Since this is taking the maximum, you should not trust this very much if you
    have specified negative potentials anywhere."""
    mm=np.zeros([mat.shape[0]//2,mat.shape[1]//2],dtype=mat.dtype)
    mm=np.maximum(mm,mat[::2,::2])
    mm=np.maximum(mm,mat[::2,1::2])
    mm=np.maximum(mm,mat[1::2,::2])
    mm=np.maximum(mm,mat[1::2,1::2])
    return mm

def upres_mat(mat):
    """A quick & dirty way to increase the resolutio of a potential matrix by a factor
    of 2.  A smarter version here would lead to faster convergence, but I have ignored that."""
    mm=np.zeros([mat.shape[0]*2,mat.shape[1]*2],dtype=mat.dtype)
    mm[::2,::2]=mat
    mm[::2,1::2]=mat
    mm[1::2,::2]=mat
    mm[1::2,1::2]=mat
    return mm

npix=512
tol=0.01
time_1=time.time()
#set boundary conditions and mask
bc=np.zeros([npix,npix])
mask=np.zeros([npix,npix],dtype='bool')
mask[0,:]=1
mask[-1,:]=1
mask[:,0]=1
mask[:,-1]=1


#add a feature in the middle of the region, held at a constant potential.
r=npix//8
x=np.arange(0,npix)
y=np.arange(0,npix)
cx=len(x)//2
cy=len(y)//2
circle_indices=(x[np.newaxis,:]-cx)**2+(y[:,np.newaxis]-cy)**2 < r**2
bc[circle_indices]=1
mask[circle_indices]=1

npass=5
#loop through the resolutions.  In this case, start with something 2**5 times coarser than the desired resolution
#solve it, and then increase the resolution by a factor of 2.  Solve the next-higher resolution problem using that
#as the starting guess.
all_masks=[None]*npass
all_bc=[None]*npass
all_rhs=[None]*npass
all_x=[None]*npass
all_masks[0]=mask
all_bc[0]=bc
for i in range(1,npass):
    all_masks[i]=deres_mat(all_masks[i-1])
    all_bc[i]=deres_mat(all_bc[i-1])


#now, make the lowest-resolution map.  First set up the RHS given the low-res boundary conditions/masks we already made.
nn=all_masks[-1].shape[0]
niter=3*nn
all_rhs[-1]=make_rhs_2d(all_bc[-1],all_masks[-1])
all_x[-1]=cg(all_rhs[-1],0*all_rhs[-1],all_masks[-1],niter)

#and now, run a loop where you increase the resolution, solve for a fixed number of iterations, then repeat until 
#you're at your desired resolution.
for i in range(npass-2,-1,-1):
    all_rhs[i]=make_rhs_2d(all_bc[i],all_masks[i])
    x0=upres_mat(all_x[i+1])
    all_x[i]=cg(all_rhs[i],x0,all_masks[i],tol)
    #plot the current-resolution potential

#Paste in the boundary conditions since we didn't solve for them in conjugate gradient.
for i in range(npass):
    all_x[i][all_masks[i]]=all_bc[i][all_masks[i]]

#Now work out our final charge distribution.
xx=all_x[0]
rho=4*xx[1:-1,1:-1];rho=rho-xx[2:,1:-1];rho=rho-xx[:-2,1:-1];rho=rho-xx[1:-1,2:];rho=rho-xx[1:-1,:-2]

time_2=time.time()
print('This entire process took ', time_2-time_1, 's')

#Finally, plot the final map with boundary conditions enforced as well as charge density
plt.figure()
plt.imshow(all_x[0])
plt.colorbar()
plt.title('Numerical V')
plt.savefig('q3_final_V.pdf')
plt.figure()
plt.imshow(rho)
plt.colorbar()
plt.title('Final charge distribution')
plt.savefig('q3_final_rho.pdf')