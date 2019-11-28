# Assignment 4

Most of the questions asked in the assignment are answered either in the comments or in the print statements shown when the code 
is ran. If you don't want to run the script, I've included screenshots of the outputs in this repo. This readme will mostly serve
to direct you to where all the answers are, or to clarify some minor things.

I'm not too sure if I ever mention this in the comments, but I assume that the noise is Gaussian and stationary, which allows
me to assume that N is diagonal and allows for a simpler way to pre-whiten, as well as find matches and snr.

## Part a)

This is all answered in comments,and smoothed power spectra are shown in plots, the filenames of which are printed when the code
is ran.

##Part b)

This is also all answered in the comments and the matches are plotted

## Part c)

idem but it is the SNRs that are plotted

## Part d)

For this part, the analytic and actual SNRs are very well lined up, but the maximum SNRs are not quite equal. This is probably
due to the fact that assuming that noise is stationary and gaussian, which is not really exact, since there are clearly some
features that repeat in the strains.

## Part e)

The answer to this is printed for each event when my code is run

## Part f)

As explained in the comments I fit a Gaussian to the peak in both SNRs to find an estimate for the arrival times and the
uncertainty on that. I then estimate the error on the position angle to be c times the time offset between the two detectors
divided by the distance between the two. These results are all printed when my code is ran.
