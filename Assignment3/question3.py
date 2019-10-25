import numpy as np
from matplotlib import pyplot as plt
import camb

#Provided function to fit
def get_spectrum(pars,lmax):
    #print('pars are ',pars)
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

#Generate step in params based on covariance matrix
def take_step_cov(covmat):
    mychol=np.linalg.cholesky(covmat)
    return np.dot(mychol,np.random.randn(mychol.shape[0]))

#Load covariance matrix and data, declare initial params and some other parameters
mycov=np.loadtxt('cov_mtx_q2.txt')
params=np.asarray([70, 2.20562123e-02, 1.21805670e-01, 5.00000000e-02, 2.08027638e-09, 9.49973440e-01])
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
x=wmap[:,0]
y=wmap[:,1]
y_err=wmap[:,2]
nstep=3000
chains=np.zeros([nstep,len(params)])
chisqvec=np.zeros(nstep)
chisq=np.sum( (y-get_spectrum(params,len(x))[2:])**2/y_err**2)
scale_fac=0.5
for i in range(nstep):
    #Try new parameters
    new_params=params+take_step_cov(mycov)*scale_fac
    new_model=get_spectrum(new_params,len(x))[2:]
    new_chisq=np.sum( (y-new_model)**2/y_err**2)
    delta_chisq=new_chisq-chisq
    prob=np.exp(-0.5*delta_chisq)
    #Decide if we switch to new params (Do not switch if tau is negative)
    accept=np.random.rand(1)<prob
    if accept and new_params[3] > 0:
        params=new_params
        model=new_model
        chisq=new_chisq
        print('Chose new fit params with chi2=', chisq, ' on iteration ', i)
    #Update chain and chi2 vector
    chains[i,:]=params
    chisqvec[i]=chisq

print('Chain Done! Run plot_chains.py to plot chains and get best params')
np.savetxt('q3_chains.txt', chains)