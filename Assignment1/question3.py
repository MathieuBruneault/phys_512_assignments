import numpy as np


sig=0.1

def fun(x):
    return 1.0/(1.0+x**2)

def fun2(x):
    return 1.0+np.exp(-0.5*x**2/(sig**2))

def simple_integrate(fun,a,b,tol):
    x=np.linspace(a,b,5)
    dx=(b-a)/4.0
    #np.median(np.diff(x))
    y=fun(x)
    neval=len(x) #let's keep track of function evaluations
    f1=(y[0]+4*y[2]+y[4])/6.0*(b-a)
    f2=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12.0*(b-a)
    myerr=np.abs(f2-f1)
    #print([a,b,f1,f2])
    if (myerr<tol):
        #return (f2)/1.0,myerr,neval
        return (16.0*f2-f1)/15.0,myerr,neval
    else:
        mid=0.5*(b+a)
        f_left,err_left,neval_left=simple_integrate(fun,a,mid,tol/2.0)
        f_right,err_right,neval_right=simple_integrate(fun,mid,b,tol/2.0)
        neval=neval+neval_left+neval_right
        f=f_left+f_right
        err=err_left+err_right
        return f,err,neval


def integrate(fun,a,b,fa,fb,fmid,tol):
    x=np.linspace(a,b,5)
    #Make y array, only needing to make 2 evaluations
    y=[fa,fun(x[1]),fmid,fun(x[3]),fb]
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
        f_left,err_left,neval_left=integrate(fun,a,mid,y[0],y[2],y[1],tol/2.0)
        f_right,err_right,neval_right=integrate(fun,mid,b,y[2],y[4],y[3],tol/2.0)
        neval=neval+neval_left+neval_right
        f=f_left+f_right
        err=err_left+err_right
        return f,err,neval

def start_integral(fun,a,b,tol,funtype,pred):
    #Integrate the function using the better integrator
    f,err,neval = integrate(fun,a,b,fun(a),fun(b),fun((a+b)/2.0),tol)

    #Dont forget to add the 3 evaluation that we made at the very beginning!
    neval += 3

    #Show results
    print('For ',funtype,' function and my integrator, f,err,neval are ', repr([f,err,neval]), ' with err ', repr(np.abs(f-pred)))




f,err,neval=simple_integrate(np.exp,-1,1,1e-3);pred=np.exp(1)-np.exp(-1)
print('For exponential function and the simple integrator, f,err,neval are ' + repr([f,err,neval])+' with err ' + repr(np.abs(f-pred)))
start_integral(np.exp,-1,1,1e-3,'exponential',pred)

f,err,neval=simple_integrate(fun,-1,1,1e-4);pred=np.arctan(1)-np.arctan(-1)
print('For fun1 function and the simple integrator, f,err,neval are ' + repr([f,err,neval])+' with err ' + repr(np.abs(f-pred)))
start_integral(fun,-1,1,1e-4,'fun1',pred)

a=-5;b=5;f,err,neval=simple_integrate(fun2,a,b,1e-4);pred=(b-a)+np.sqrt(2*np.pi)*sig
print('For gaussian function and the simple integrator, f,err,neval are ' + repr([f,err,neval])+' with err ' + repr(np.abs(f-pred)))
start_integral(fun2,a,b,1e-4,'gaussian',pred)
