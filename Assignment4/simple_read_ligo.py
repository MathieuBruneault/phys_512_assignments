import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import matplotlib.mlab as mlab

def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    th=template[0]
    tl=template[1]
    return th,tl
def read_file(filename):
    dataFile=h5py.File(filename,'r')
    dqInfo = dataFile['quality']['simple']
    qmask=dqInfo['DQmask'][...]
    meta=dataFile['meta']
    gpsStart=meta['GPSstart'][()]
    utc=meta['UTCstart'][()]
    duration=meta['Duration'][()]
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)
    dataFile.close()
    return strain,dt,utc
def get_fft_vec(n):
    vec=np.arange(n)
    vec[vec>n/2]=vec[vec>n/2]-n
    return vec
def smooth_map(map,npix,smooth=True):
    nx=len(map)
    xind=get_fft_vec(nx)

    #make a 1-d gaussians of the correct length
    sig=npix/np.sqrt(8*np.log(2))
    xvec=np.exp(-0.5*xind**2/sig**2)/np.sqrt(2*np.pi*sig**2)  
    xvec=xvec/xvec.sum()

    #if we didn't mess up, the kernel FT should be strictly real
    kernel=np.real(np.fft.rfft(xvec))
    
    #get the map Fourier transform
    mapft=np.fft.rfft(map)
    #Multiply it by kerel to smooth
    mapft=mapft*kernel

    return mapft

def noise_model(fname):
    strain,dt,utc=read_file(fname)
    x=np.arange(len(strain))
    x=x-1.0*x.mean()
    window=0.5*(1+np.cos(x*np.pi/np.max(x)))
    windowed_strain=window*strain
    normfac=np.sqrt(np.mean(window**2))
    power_spectrum=(np.abs(smooth_map(windowed_strain,1)/normfac)**2)
    power_spectrum=power_spectrum[5000:]
    power_spectrum=power_spectrum[:-15000]
    plt.figure()
    plt.plot(power_spectrum)
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    return power_spectrum

def find_wave(fname_H, fname_L, template_name):
    strain_H,dt,utc=read_file(fname_H)
    strain_L,dt,utc=read_file(fname_L)
    template_H,template_L=read_template(template_name)
    noise_H=np.fft.fftshift(np.fft.irfft(noise_model(fname_H)))
    noise_L=np.fft.fftshift(np.fft.irfft(noise_model(fname_L)))
    whitened_strain_H=((noise_H)**(-1/2))*strain_H
    whitened_strain_L=((noise_L)**(-1/2))*strain_L
    whitened_template_H=((noise_H)**(-1/2))*template_H
    whitened_template_L=(noise_L)**(-1/2)*template_L
    plt.figure()
    plt.plot(whitened_template_H)
    plt.show()



fname_H='LOSC_Event_tutorial/H-H1_LOSC_4_V2-1126259446-32.hdf5'
fname_L='LOSC_Event_tutorial/L-L1_LOSC_4_V2-1126259446-32.hdf5'
template_name='LOSC_Event_tutorial/GW150914_4_template.hdf5'
find_wave(fname_H, fname_L, template_name)





































































































































































































































































































































































































































































































































































































































































































































































































































































