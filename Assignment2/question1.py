import numpy as np
from matplotlib import pyplot as plt

def fun1(x):
    #Log function
    return np.log2(x)

def cheby_fit(x,y,ord):
    x_new = np.zeros(len(x))

    #Map the given interval onto -1,1
    for i in range(len(x)):
        x_new[i] =  -1 + ((x[i] - x[0]) * 2 / (x[-1] - x[0]))

    #Generate the A matrix
    mat=np.zeros([len(x),ord+1])
    mat[:,0]=1.0
    if ord>0:
        mat[:,1]=x_new
    if ord>1:
        for i in range(1,ord):
            mat[:,i+1]=2*x_new*mat[:,i]-mat[:,i-1]

    #Compute the fit parameters
    lhs=np.dot(mat.transpose(),mat)
    rhs=np.dot(mat.transpose(),y)
    fitp=np.dot(np.linalg.inv(lhs),rhs)
    return fitp, mat


def get_order_and_params(x, fun, tol):
    y_true=fun(x)

    #Choose some big order and cheby fit to it
    norder=30
    fitp,mat=cheby_fit(x,y_true,norder)
    order=0

    #Truncate the fit at different order until the error is sufficiently small
    for i in range(len(fitp)):
        #err can be estimated by the last term in our truncated series
        err=np.abs(fitp[i])
        print('Chebyshev fit truncated at order ', i, ' yielded err: ', err)
        if err < tol:
            order=i
            break
    
    #Show the number of necessary terms and return the fit to the order found
    print('For this fit, ', order+1, ' terms are necessary')
    y=np.dot(mat[:,:order],fitp[:order])
    return order, y

def polyfit_and_compare(x, tol, filename):
    y_true=fun1(x)

    #Get the necessary truncation order to get our desired error
    order,y_cheb=get_order_and_params(x,fun1,tol)

    #Fit to the same order with normal polynomials
    polyparams=np.polyfit(x,y_true,order-1)
    y_poly=np.polyval(polyparams,x)

    #Compute residuals for both fits and plot away!
    res_cheb=y_cheb-y_true
    res_poly=y_poly-y_true
    plt.figure()
    plt.plot(x,res_cheb, label='Chebychev residuals')
    plt.plot(x,res_poly, label='Polynomial residuals')
    plt.legend()
    plt.title('Residuals for both fits')
    plt.xlabel('x')
    plt.ylabel('Residuals')
    plt.savefig(filename)

    #Find the max error and the rms error for each fit and print them
    cheb_max=np.amax(np.abs(res_cheb))
    cheb_rms=np.sqrt(np.mean(res_cheb**2))
    poly_max=np.amax(np.abs(res_poly))
    poly_rms=np.sqrt(np.mean(res_poly**2))
    print('The max error for Chebyshev is: ', cheb_max, ' and the rms is: ', cheb_rms)
    print('The max error for Polyfit is: ', poly_max, ' and the rms is: ', poly_rms)

#Fit and compare cheby vs. polyfit for the ranges 0.5,1 and 0.5,3
print('Starting fit from 0.5 to 1')
x=np.linspace(0.5,1,1000)
polyfit_and_compare(x,1e-6,'question1_fit_residuals_half_to_1.pdf')
print('Starting fit from 0.5 to 3')
x=np.linspace(0.5,3,1000)
polyfit_and_compare(x,1e-6,'question1_fit_residuals_half_to_3.pdf')

#NOTE: THE CODE IS MADE SUCH THAT IT ANSWER BOTH A AND B