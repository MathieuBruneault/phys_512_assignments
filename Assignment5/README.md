# Assignment 5

##Question 1

The numerical V and charge density are repectively shown in q1_final_V.pdf and q1_final_rho.pdf. The charge looks all bundled
up on the surface of the wire, as should be expected. The charge per unit length was found to be 0.00016996735801284984, which
is found by summing up all the charge and dividing by the circumference of our cylinder, assuming that all the charge is indeed
on its surface. The analytic solution is presented in q1_true_V.pdf. It might not be super exact because we are kind of ignoring
boundary conditions here, nonetheless the true and numerical solutions look pretty similar.

##Question 2

In question 1, finding V for a tolerance of 0.01 took 10488 steps and 45.534708738327026 s, where it only took 85 steps and
0.682905912399292 s to run using the conjugate method. Using this method indeed saved a lot of time.
The numerical V and charge density are repectively shown in q2_final_V.pdf and q2_final_rho.pdf.
They look a bit different, but that might be because tolerance is a bit low, given that we needed a resonnable tolerance 
for #1 to not take ages.

##Question 3

As mentionned, in qusetion 2 we got a numerical V in 0.682905912399292 s. In this part, using the changes of resolution, we got
V up to the same tolerance in only 0.16097283363342285 s, so it indeed got faster again.
The numerical V and charge density are repectively shown in q3_final_V.pdf and q3_final_rho.pdf.

##Question 4

The numerical V (with E) and charge density with a bump are repectively shown in q4_final_V_and_E.pdf and q4_final_rho.pdf.
As it is not very obvious what is wrong just from looking at the E field in 4_final_V_and_E.pdf, I've also included a plot
of the strength of the E field as a function of y for x fixed in the middle of the box, corresponding to the middle of the
wire and the bump. We can see that for E just before getting to the cylinder (where E is 0 as expected), there is a weird
spike which does not happen on the other side, ehwhere there is no bump. So the bump indeed disrupts E.

##Question 5

For this question we are told to plot T along a line from the center of a box tothe center on the other side. This is just 
the exact same as solving the heat equation in 1D so that's what I did. Also, we aren't told what the boundary conditions
are on the right wall, hence I assume either Dirichlet or Neumann is fine so I chose Neumann (for no particular reason).
I you run the script q5.py it'll plot the temperature distribution a a bunch of different times and show all of the. It 
also saves a figure of what the distribution looked like at 3 different times. These are shown in the q5_temp_dist files.
I chose k=10^-4, which yields nice solutions. If k is an order of magnitude hgher or more, we seem to be in a different
regime, and the solutions aren't as nice
