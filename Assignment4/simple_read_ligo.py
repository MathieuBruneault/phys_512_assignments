import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
import matplotlib.mlab as mlab
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
import os
plt.rcParams.update({'figure.max_open_warning': 0})

path_to_data='LOSC_Event_tutorial/'

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

def noise_model(window, strain):
    #Window the strains using a cosine function going to 0 and its norm factor
    windowed_strain=window*strain
    normfac=np.sqrt(np.mean(window**2))

    #Estimate the power spectrum as the square of the FT of Strains
    power_spectrum=np.abs(np.fft.rfft(windowed_strain)/normfac)**2

    #Apply a gaussian filter to smooth the power spectrum
    power_spectrum=gaussian_filter(power_spectrum,1)

    return power_spectrum
def gaussian(x,A,sig,mu):
    return A*np.exp(-(x-mu)**2/sig**2)

def find_wave(event_name, fname_H, fname_L, template_name):
    print('------------------Processing event ' + event_name + ' -----------------------')
    print('Part a)')
    #Get the strains, template, corresponding times and frequencies for both detectors
    strain_H,dt,utc_H=read_file(fname_H)
    strain_L,dt,utc_L=read_file(fname_L)
    time=np.linspace(0,len(strain_H)/4096,len(strain_H))
    freq=np.fft.rfftfreq(len(strain_H),dt)
    template_H,template_L=read_template(template_name)

    #Calculate a window function to window our data
    x=np.arange(len(strain_H))
    x=x-1.0*x.mean()
    window=0.5*(1+np.cos(x*np.pi/np.max(x)))

    #Don't forget the norm factor
    normfac=np.sqrt(np.mean(window**2))

    #Get the power spectrum, which will be our noise model
    noise_H=noise_model(window, strain_H)
    noise_L=noise_model(window, strain_L)

    #Plot the power spectra in log log scale
    filename_H=event_name + '/PS_H.pdf'
    filename_L=event_name + '/PS_L.pdf'
    plt.figure()
    plt.loglog(freq,noise_H)
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Power spectrum')
    plt.savefig(filename_H)
    plt.figure()
    plt.xlabel('Freq (Hz)')
    plt.ylabel('Power spectrum')
    plt.loglog(freq,noise_L)
    plt.savefig(filename_L)
    print('The power spectra (corresponding to my noise models in Frequency space) are plotted in the files ' + filename_H + ' and ' + filename_L)

    print('Part b)')
    #Pre-whiten the template and the strains (not forgetting to window beforehand) to get A and d
    #in a basis where N is the identity matrix
    whitened_A_H=np.fft.rfft(window*template_H)/(np.sqrt(noise_H)*normfac)
    whitened_A_L=np.fft.rfft(window*template_L)/(np.sqrt(noise_L)*normfac)
    whitened_d_H=np.fft.rfft(window*strain_H)/(np.sqrt(noise_H)*normfac)
    whitened_d_L=np.fft.rfft(window*strain_L)/(np.sqrt(noise_L)*normfac)

    #Apply matched filter using the noise model found previously and convolving our whitened A and d
    m_H=np.fft.fftshift(np.fft.irfft(np.conj(whitened_A_H)*whitened_d_H))
    m_L=np.fft.fftshift(np.fft.irfft(np.conj(whitened_A_L)*whitened_d_L))

    #Plot the matched filters
    filename_H=event_name + '/Match_H.pdf'
    filename_L=event_name + '/Match_L.pdf'
    plt.figure()
    plt.plot(time,m_H)
    plt.xlabel('Time (s)')
    plt.ylabel('Match')
    plt.savefig(filename_H)
    plt.figure()
    plt.xlabel('Time (s)')
    plt.ylabel('Match')
    plt.plot(time,m_L)
    plt.savefig(filename_L)
    print('The matches are plotted in time domain in the files ' + filename_H + ' and ' + filename_L)

    print('Part c)')
    #Find the SNR for both detectors, as well as their combined SNR. The noise estimate is the power spectrum in time domain
    SNR_H=np.abs(m_H*np.fft.fftshift(np.fft.irfft(np.sqrt(np.conj(whitened_A_H)*whitened_A_H))))
    SNR_L=np.abs(m_L*np.fft.fftshift(np.fft.irfft(np.sqrt(np.conj(whitened_A_L)*whitened_A_L))))
    SNR_Comb=np.sqrt(SNR_H**2+SNR_L**2)

    #Plot the SNRs
    filename=event_name + '/SNRs.pdf'
    fig, axs = plt.subplots(3)
    axs[0].plot(time,SNR_H)
    axs[1].plot(time,SNR_L)
    axs[2].plot(time,SNR_Comb)
    axs[0].set(ylabel='SNR Hanford')
    axs[1].set(ylabel='SNR Livingston')
    axs[2].set(ylabel='SNR Combined')
    axs[2].set(xlabel='Time(s)')    
    plt.savefig(filename)

    #Print the max SNRs found
    print('The max SNR found for Hanford is ' + str(np.amax(SNR_H)))
    print('The max SNR found for Livingston is ' + str(np.amax(SNR_L)))
    print('The max combined SNR found is ' + str(np.amax(SNR_Comb)))
    print('The SNR plots are presented in ' + filename)

    print('Part d)')
    #Find the analytic by dividing the teplate by our noise model (in our case, that information is already
    #contained in our pre-whitened A, just need it back in time domain)
    Anal_SNR_H=np.abs(np.fft.irfft(whitened_A_H))
    Anal_SNR_L=np.abs(np.fft.irfft(whitened_A_L))
    Anal_SNR_Comb=np.sqrt(Anal_SNR_H**2+Anal_SNR_L**2)

    #Plot the SNRs
    filename=event_name + '/Analytic_SNRs.pdf'
    fig, axs = plt.subplots(3)
    axs[0].plot(time,Anal_SNR_H)
    axs[1].plot(time,Anal_SNR_L)
    axs[2].plot(time,Anal_SNR_Comb)
    axs[0].set(ylabel='Analytic SNR Hanford')
    axs[1].set(ylabel='Analytinc SNR Livingston')
    axs[2].set(ylabel='Analytic SNR Combined')
    axs[2].set(xlabel='Time(s)')
    plt.savefig(filename)

    #Print the max SNRs found
    print('The max SNR found for Hanford is ' + str(np.amax(Anal_SNR_H)))
    print('The max SNR found for Livingston is ' + str(np.amax(Anal_SNR_L)))
    print('The max combined SNR found is ' + str(np.amax(Anal_SNR_Comb)))
    print('The SNR plots are presented in ' + filename)

    print('Part e)')
    #Find the cummulative sum of power spectra
    PS_sum_H=np.cumsum(np.abs(whitened_A_H**2))
    PS_sum_L=np.cumsum(np.abs(whitened_A_L)**2)

    #Find the frequency at which this cummulative sum is closest to half of its maximum
    half_power_freq_H=freq[np.argmin(np.abs(PS_sum_H-(np.amax(PS_sum_H)/2)))]
    half_power_freq_L=freq[np.argmin(np.abs(PS_sum_L-(np.amax(PS_sum_L)/2)))]
    print('For this event, the half power frequency is ' + str(half_power_freq_H) + ' Hz for Hanford and ' +
    str(half_power_freq_L) + ' Hz for Livingston')

    print('Part f)')
    #To get an estimate on arrival time and its error we fit a gaussian to each SNR vs time plot
    #and fit a Gaussian. The mean of the Gaussian the arrival time and sigma is the error.
    A_H=np.amax(SNR_H)
    mu_H=time[np.argmax(SNR_H)]
    sig_H=0.001
    params_H,cov=curve_fit(gaussian, time[np.argmax(SNR_H)-5:np.argmax(SNR_H)+5], SNR_H[np.argmax(SNR_H)-5:np.argmax(SNR_H)+5], p0=[A_H,sig_H,mu_H])
    A_L=np.amax(SNR_L)
    mu_L=time[np.argmax(SNR_L)]
    sig_L=0.001
    params_L,cov=curve_fit(gaussian, time[np.argmax(SNR_L)-5:np.argmax(SNR_L)+5], SNR_L[np.argmax(SNR_L)-5:np.argmax(SNR_L)+5], p0=[A_L,sig_L,mu_L])

    #Print the arrival times and the errors on them
    print('The arrival time found for Hanford is ' + str(params_H[2]) + ' +/- ' + str(params_H[1]) + ' seconds after the initial time of ' + str(utc_H))
    print('The arrival time found for Livingston is ' + str(params_L[2]) + ' +/- ' + str(params_L[1]) + ' seconds after the initial time of ' + str(utc_L))
    
    

#Create a dict with all file names
events=[
    {
        "name":"GW150914",
        "fn_H1"       : "H-H1_LOSC_4_V2-1126259446-32.hdf5",
        "fn_L1"       : "L-L1_LOSC_4_V2-1126259446-32.hdf5",
        "fn_template" : "GW150914_4_template.hdf5"
    },
    {
        "name":"LVT151012",
        "fn_H1"       : "H-H1_LOSC_4_V2-1128678884-32.hdf5",
        "fn_L1"       : "L-L1_LOSC_4_V2-1128678884-32.hdf5",
        "fn_template" : "LVT151012_4_template.hdf5"
    },
    {
        "name":"GW151226",
        "fn_H1"       : "H-H1_LOSC_4_V2-1135136334-32.hdf5",
        "fn_L1"       : "L-L1_LOSC_4_V2-1135136334-32.hdf5",
        "fn_template" : "GW151226_4_template.hdf5"
    },
    {
        "name":"GW170104",
        "fn_H1"       : "H-H1_LOSC_4_V1-1167559920-32.hdf5",
        "fn_L1"       : "L-L1_LOSC_4_V1-1167559920-32.hdf5",
        "fn_template" : "GW170104_4_template.hdf5"
    }
    ]

#Loop over all events
for event in events:
    directory = './' + event["name"]
    if not os.path.exists(directory):
        os.makedirs(directory)
    find_wave(event["name"], path_to_data + event["fn_H1"], path_to_data + event["fn_L1"], path_to_data + event["fn_template"])




































































































































































































































































































































































































































































































































































































































































































































































































































































