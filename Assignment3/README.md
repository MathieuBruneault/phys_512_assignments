# Assignment 3

## Question 1

With the initial parameters given I got Chi^2=1588.436665

## Question 2

With the parameters given and a fixed tau=0.05, Newton's method converged succesfully, hence there was
no need to use Levenberg-Marquadt algorithm. Here are the best fit parameters and errors I got:

* H0 = 63 +/- 3
* wb h^2 = 0.0221 +/- 0.0008
* wc h^2 = 0.122 +/- 0.007
* tau = 0.05 +/- 0.15
* As = (2.1 +/- 0.6) x 10^-9
* slope = 0.95 +/- 0.03

These errors were computed using a covariance matrix, so the huge uncertainty on tau indicates that it
is not very correlated with the other parameters. This means that if I were to run this with float tau,
not much information would be gained on the other parameters by fitting, so the errors would not change much.

The derivatives I use are computed numerically. My code prints out values of Chi^2 on each iteration, and
since it keeps getting better before finding a minimum at around Chi^2=1245, we can assume my estimates
for the dervatives are pretty accurate, otherwise Newton's method would just not work.

## Question 3

Here are the parameters and their errors I get with mcmc:

* H0 = 70 +/- 3
* wb h^2 = 0.0227 +/- 0.0007
* wc h^2 = 0.112 +/- 0.005
* tau = 0.09 +/- 0.06
* As = (2.2 +/- 0.2) x 10^-9
* slope = 0.98 +/- 0.02

Looking at the plots in q3_chains.pdf we see that the chains mostly look like noise by the end of the
3000 iterations, indicating that they have converged by then.

## Question 4

I added Planck's tau value as a prior by rejecting all parameters where tau wasn't within 3 sigma of 
that value, so outside of their 99% certainty region. Here are the full results of that chain:

* H0 = 69 +/- 3
* wb h^2 = 0.0222 +/- 0.0006
* wc h^2 = 0.114 +/- 0.007
* tau = 0.054 +/- 0.013
* As = (2.06 +/- 0.07) x 10^-9
* slope = 0.965 +/- 0.014

These are all consistent within 1 sigma with the values obtained in question 3 and the errors are
pretty similar, except the error on tau, which got smaller.

By importance sampling chains from question 3, I got the following parameters

* H0 = 69 +/- 3
* wb h^2 = 0.0221 +/- 0.0006
* wc h^2 = 0.112 +/- 0.007
* tau = 0.054 +/- 0.003
* As = (2.04 +/- 0.07) x 10^-9
* slope = 0.963 +/- 0.014

This yields results consistent within 1 sigma of both chains in question 3 and question 4
