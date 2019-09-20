import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

#This is taken from solutions to Griffiths 2.7 and needs to be integrated for x from -1 to 1.
#All constants (including R) are set to 1
def electric_field(x,z):
    return (z-x)/((1+z**2-2*z*x)**(3/2))

def my_integrate(fun,a,b,fa,fb,fmid,tol,z):
    x=np.linspace(a,b,5)
    #Make y array, only needing to make 2 evaluations
    y=[fa,fun(x[1],z),fmid,fun(x[3],z),fb]
    neval=2
    #Calculate f1, f2 and err as in the simple integrator
    f1=(y[0]+4*y[2]+y[4])/6.0*(b-a)
    f2=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12.0*(b-a)
    myerr=np.abs(f2-f1)
    if (myerr<tol):
        #Return relevant values if success
        return (16.0*f2-f1)/15.0,myerr,neval
    else:
        mid=0.5*(b+a)

        #Call the integrator again, providing boundary and middle y values (which were already evaluated in y)
        #to avoid calling fun unnecessarily.
        f_left,err_left,neval_left=my_integrate(fun,a,mid,y[0],y[2],y[1],tol/2.0,z)
        f_right,err_right,neval_right=my_integrate(fun,mid,b,y[2],y[4],y[3],tol/2.0,z)
        neval=neval+neval_left+neval_right
        f=f_left+f_right
        err=err_left+err_right
        return [f,err,neval]

#This includes z points from 0 to 2R and one of the points is z = R
z_range = np.linspace(0,2,101)
quad_int = []
my_int = []

for z in z_range:
    quad_int = np.append(quad_int,integrate.quad(electric_field,-1,1,args=(z))[0])
    my_int = np.append(my_int,my_integrate(electric_field,-1,1,electric_field(-1,z),electric_field(1,z),electric_field(0,z),1e-3,z)[0])

plt.figure()
plt.plot(z_range, quad_int)
plt.savefig('quad_integral.pdf')

plt.figure()
plt.plot(z_range, my_int)
plt.savefig('my_integral.pdf')

#This code will crash as it is, because my integrator cannot handle the singularity at z=1 (corresponding to z=R)
#Quad does not care about the singularity and produces the plot quad_integral.pdf even if one of my z is R.
#To see quad produce results, you can comment out lines 47 to 49 and line 41.
#If the value z=R is not included in the range (for instance by replacing the 101 at line 35 by 150), then
#my integrator does work and produces the plot my_integral.pdf