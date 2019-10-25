import numpy as np
import camb
from matplotlib import pyplot as plt
import time

#Provided function using Camb for our model
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


plt.figure()

#Declare initial parameters
pars=np.asarray([65,0.02,0.1,0.05,2e-9,0.96])
#Get data and its length
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
lmax=len(wmap[:,0])

plt.clf()
#Plot our data
plt.plot(wmap[:,0],wmap[:,1],'.')

#Execute fit with basic parameters, get chi2 and plot fit
cmb=get_spectrum(pars, lmax)[2:]
chi2=np.sum(((wmap[:,1]-cmb)**2)/wmap[:,2]**2)
print('Chi^2 obtained with initial parameters is',chi2)
plt.plot(cmb)

plt.show()