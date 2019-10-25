import numpy as np
from matplotlib import pyplot as plt
import camb

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

def take_step_cov(covmat):
    mychol=np.linalg.cholesky(covmat)
    return np.dot(mychol,np.random.randn(mychol.shape[0]))

def mcmc(wmap):
    #This is mostly the same as in question3.py, with a prior on tau
    mycov=np.loadtxt('cov_mtx_q2.txt')
    params=np.asarray([70, 2.20562123e-02, 1.21805670e-01, 5.00000000e-02, 2.08027638e-09, 9.49973440e-01])
    wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
    x=wmap[:,0]
    y=wmap[:,1]
    y_err=wmap[:,2]
    nstep=3000
    chains=np.zeros([nstep,len(params)]) #keep track of where the chain went
    chisqvec=np.zeros(nstep)
    chisq=np.sum( (y-get_spectrum(params,len(x))[2:])**2/y_err**2)
    scale_fac=0.5
    for i in range(nstep):
        new_params=params+take_step_cov(mycov)*scale_fac
        #Only difference from 3: Only accept new params if tau is it's in
        #in the 99% certaity region of Plank's tau value
        if new_params[3]<0.0325 or new_params[3]>0.0763:
            chains[i,:]=params
            chisqvec[i]=chisq
            continue
        new_model=get_spectrum(new_params,len(x))[2:]
        new_chisq=np.sum( (y-new_model)**2/y_err**2)
        delta_chisq=new_chisq-chisq
        prob=np.exp(-0.5*delta_chisq)
        accept=np.random.rand(1)<prob
        if accept and new_params[3] > 0:
            params=new_params
            model=new_model
            chisq=new_chisq
            print('Chose new fit params with chi2=', chisq, ' on iteration ', i)
        chains[i,:]=params
        chisqvec[i]=chisq
    print('Chain Done! Run plot_chains.py to plot chains and get best params')
    return chains, chisqvec

#Load in the chains from Q3
chains_q3=np.loadtxt('q4_chains.txt')
nsteps=len(chains_q3[:,0])
#Calculate weights, weighing using Planck's value for tau
weights=1/(chains_q3[int(nsteps/2):,:]-0.0544)**2
params=np.average(chains_q3[int(nsteps/2):,:],weights=weights,axis=0)
#Estimate error using weighted variance
variance=np.average((chains_q3[int(nsteps/2):,:]-params)**2,weights=weights,axis=0)
err=np.sqrt(variance)
print('Best fit params are ', params)
print('Errors are: ', err)

#Run mcmc with prior on tau
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
chains, chisqvec = mcmc(wmap)
np.savetxt('q4_chains.txt', chains)