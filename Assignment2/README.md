# Assignment 2

This readme will answer the questions asked in assignment 2

## Question 1 a)

Using a truncated Chebyshev fit (truncated from order 30), it is necessary to go up to order 8 to get an accuracy
better than 10^-6, so I need 9 terms to get that accuracy. A polynomial fit of the same order is carried out, and
the residuals are plotted in the file question1_fit_residuals_half_to_1.pdf. The max errors compare as follows for
both fits:

* Max Error = 3.1969782843255246e-07 for truncated Chebyshev fit
* Max Error = 7.888699600222537e-07 for a polynomial fit of the same order'

And the rms compare as follows:

* rms error = 1.9196349459354577e-07 for truncted Chebyshev fit'
* rms error = 1.685251825365508e-07 for a polynomial fit of the same order

## Question 1 b)

The code from question 1 a) was made to be general enough to also answer this question. The residuals for a fit
of log base 2 from 0.5 to 3 are plotted in the file question1_fit_residuals_half_to_3.pdf.

## Question 2 a)

The model I will use is a decaying exponential A*e^(-(t-t0)/tau)+C. It would be possible to include the A in the
t0 and tau constants, but the in the model I use, all the parameters have nice physical interpretations: A is the
amplitude, t0 is the starting time of the flare, tau is the time constant and C is the offset. This model is 
non-linear, since it includes an exponential. This model looks like an exponential decay. The initial parameters are:

* A=0.25, since the amplitude of the exponential seems to approximately be 0.25
* t0=1706.52314133376, the time where flux is maximal, so when the flare starts
* tau=0.01, found by trying different values and plotting the guess against data  
* C=1, seems to be the approximate background

The data and my initial guess, zoomed in around t=1706, are plotted in the file question2_initial_fit.pdf

## Question 2 b)

Newton's method yields the following fit parameters:

* t0 = 1706.52118820876
* tau = 0.015332728236290921
* A = 0.296875
* C = 1.0020997204120747

The best fit is plotted along the initial fit and the data in the file question2_initial_and_final_fits.pdf

## Question 2 c)

I model my errors by generating 50 sets of fake data and fitting to it. Taking my best fit parameters, I generate
fake data by adding random gaussian noise to my best fit, taking the standard deviation of this noise to be the 
rms error on between my best fit and the data. This yields 50 sets of 'simulated' parameters, and I take the 
error on my best fit parameters to be the standard deviation of these 'simulated' parameters. This yields the 
following errors:

* t0 err = 0.00193516
* tau err = 0.00049156
* A err = 0.03759786
* C err = 0.00071523

This seems reasonable, the errors on t0 and C are small compared to their associated parameters because they are 
more well constrained, while tau and A aren't as well constrained, so they have slightly bigger % error.

## Question 2 d)

In the way I model my errors, I assume gaussian noise. Looking at the full span of data, there is clearly some
correlated noise, which has not been taken into consideration in this model. Hence, I would not trust these
errors too much.
