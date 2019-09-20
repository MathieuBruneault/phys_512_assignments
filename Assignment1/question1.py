import numpy as np
import matplotlib.pyplot as plt

eps=1e-16
#Define the functions we will need to use and their derivatives

def fun1(x):
    return np.exp(x)

def fun2(x):
    return np.exp(0.01*x)

def first_deriv_fun1(x):
    return np.exp(x)

def first_deriv_fun2(x):
    return 0.01*np.exp(0.01*x)

#We will test whether our estimate for delta really is ideal for x=0 by finding
#the error in the approximation for a range of deltas
def find_error(fun, fun_deriv, delta):
    approx_deriv = (8*(fun(delta)-fun(-delta))-(fun(2*delta)-fun(-2*delta)))/(12*delta)
    actual_deriv = fun_deriv(0)
    err = np.abs(approx_deriv-actual_deriv)
    return err

#Find the ideal delta for both functions
fun1_del = ((45*fun1(0)*eps)/(8*fun1(0)))**(1/5)
fun2_del = ((45*fun2(0)*eps)/(8*(0.01**5)*fun2(0)))**(1/5)

#Find the error for a range of deltas for both functions
fun1_x = np.linspace(0.001*fun1_del, 4*fun1_del, 200)
fun1_err = []
for delta in fun1_x:
    fun1_err = np.append(fun1_err, find_error(fun1, first_deriv_fun1, delta))

fun2_x = np.linspace(0.001*fun2_del, 4*fun2_del, 200)
fun2_err = []
for delta in fun2_x:
    fun2_err = np.append(fun2_err, find_error(fun2, first_deriv_fun2, delta))

plt.figure()
plt.plot(fun1_x, fun1_err)
plt.plot([fun1_del, fun1_del], [0,np.amax(fun1_err)], 'k', label='Ideal Delta')
plt.ylabel('Error')
plt.xlabel('Delta')
plt.title('Error in the derivative approximaton vs delta for exp(x) at x=0')
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig('error_analysis_fun1.pdf')

plt.figure()
plt.plot(fun2_x, fun2_err)
plt.plot([fun2_del, fun2_del], [0,np.amax(fun2_err)], 'k', label='Ideal Delta')
plt.ylabel('Error')
plt.xlabel('Delta')
plt.title('Error in the derivative approximaton vs delta for exp(0.01x) at x=0')
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig('error_analysis_fun2.pdf')

#From the two plots produced, we can see that the ideal edlta more or less aligns with the minimum in error