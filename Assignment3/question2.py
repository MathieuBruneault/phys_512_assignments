import numpy as np
from matplotlib import pyplot as plt
import camb

#Provided function to fit
def get_spectrum(pars,lmax, fixTau):
    if fixTau:
        #Reinsert tau=0.05 in params
        pars=np.insert(pars,3,0.05)
    H0=pars[0]
    ombh2=pars[1]
    omch2=pars[2]
    tau=pars[3]
    As=pars[4]
    ns=pars[5]
    pars=camb.CAMBparams()
    pars.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    pars.InitPower.set_params(As=As,ns=ns,r=0)
    pars.set_for_lmax(lmax,lens_potential_accuracy=0)
    results=camb.get_results(pars)
    powers=results.get_cmb_power_spectra(pars,CMB_unit='muK')
    cmb=powers['total']
    tt=cmb[:,0]
    return tt

def get_grad(pars,lmax,fixTau):
    #Get the gradient by numerically computing it
    grad=np.zeros([lmax,len(pars)])
    #Loop through parameters to compute the numerical derivative for each
    for i, par in enumerate(pars):
        #Get an appropriate dx depending on the parameter
        dx=par*0.01
        #Compute the numerical derivative using f(x-dx) and f(x+dx)
        pars_below=pars.copy()
        pars_below[i]-=dx
        f_below=get_spectrum(pars_below,lmax,fixTau)[2:]
        pars_above=pars.copy()
        pars_above[i]+=dx
        f_above=get_spectrum(pars_above,lmax,fixTau)[2:]
        deriv=(f_above-f_below)/(2*dx)
        grad[:,i]=deriv
    return grad

def newton(wmap,p0,fixTau):
    p=p0.copy()
    if fixTau:
        #Remove tau from params to not fit for it
        p=np.delete(p,3)
    #Run Newton's Method as usual and print Chi2 at each iteration
    for j in range(5):
        pred=get_spectrum(p,len(wmap[:,0]),fixTau)[2:]
        grad=get_grad(p,len(wmap[:,0]),fixTau)
        r=np.array(wmap[:,1]-pred)
        chi2=(r**2/wmap[:,2]**2).sum()
        print('On iteration ', j, ', Chi^2 was ', chi2)
        r=np.matrix(r).transpose()
        grad=np.matrix(grad)
        lhs=grad.transpose()*grad
        rhs=grad.transpose()*r
        dp=np.linalg.inv(lhs)*(rhs)
        for jj in range(len(p)):
            p[jj]=p[jj]+float(dp[jj])
    return p

def find_error(wmap,p):
    #Compute covariance matrix using best fit and get errors on params from it
    grad=get_grad(p,len(wmap[:,0]),False)
    n=np.diag(1/(wmap[:,2])**2)
    cov=np.linalg.inv(np.dot(np.dot(grad.transpose(),n),grad))
    err=np.sqrt(np.diag(cov))
    return err, cov

#Declare initial parameters and get data
p0=np.asarray([65,0.02,0.1,0.05,2e-9,0.96])
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
fixTau=True
#Get best fit and best fit params using Newton's method
best_p=newton(wmap,p0,fixTau)
best_fit=get_spectrum(best_p,len(wmap[:,0]),fixTau)[2:]

if fixTau:
    #Reinsert tau in parameters if it was not fit for
    best_p=np.insert(best_p,3,0.05)

#Find errors on params and covariance matrix
err,cov=find_error(wmap,best_p)
print('Best fit params are ', best_p)
print('Errors are ',err)

#Plot data and best fit
plt.figure()
plt.plot(wmap[:,0],wmap[:,1],'.')
plt.plot(best_fit)
plt.show()

#Save covariance matrix for q3 and q4
np.savetxt('cov_mtx_q2.txt', cov)