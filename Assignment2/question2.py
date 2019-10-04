import numpy as np
import matplotlib.pyplot as plt

def exp_model(t,p):
    #The model chosen to fit our data
    return p[2]*np.exp(-(t-p[0])/p[1])+p[3]

def grad_model(p,t):
    grad=np.zeros([len(t),len(p)])

    #Differentiate w.r.t each parameter, populate and return grad
    grad[:,0]=p[2]/p[1]*np.exp(-(t-p[0])/p[1])
    grad[:,1]=np.multiply(p[2]*(t-p[0])/(p[1]**2),np.exp(-(t-p[0])/p[1]))
    grad[:,2]=np.exp(-(t-p[0])/p[1])
    grad[:,3]=1.0
    return grad

def newton(t,y,p0):
    p=p0.copy()
    #Apply Newton's method using the gradient and return the best-fit parameters
    for j in range(5):
        pred=exp_model(t,p)
        grad=grad_model(p,t)
        r=np.array(y-pred)
        err=(r**2).sum()
        r=np.matrix(r).transpose()
        grad=np.matrix(grad)
        lhs=grad.transpose()*grad
        rhs=grad.transpose()*r
        dp=np.linalg.inv(lhs)*(rhs)
        for jj in range(len(p)):
            p[jj]=p[jj]+float(dp[jj])
        if err<1 and j!=0:
            break
    return p

def findError(t,y,p):
    #Find the rms 
    y_model=exp_model(t,p)
    rms=np.sqrt(np.mean((y-y_model)**2))
    simulated_params=[]
    #Simulate our data with random normal noise 50 times, finding new parameters each time
    #corresponding (hopefully) to redoing the experiment and re fitting.
    for i in range(50):
        noise=np.random.normal(0,rms,len(t))
        y_simulated=y_model+noise
        p_simulated=newton(t,y_simulated,p)
        simulated_params.append(p_simulated)
    #Compute our estimate for the error by taking the std of our simulated parameters
    simulated_params=np.array(simulated_params)
    errors=np.std(simulated_params,axis=0)
    return errors


data=np.loadtxt('229614158_PDCSAP_SC6.txt',delimiter=',')
#We want to fit to the highest peak so we start at the point of highest flux
idx_max=np.argmax(data[:,1])
t=data[:,0][idx_max:idx_max+50]
y=data[:,1][idx_max:idx_max+50]
#Set the initial fit parameters, explained in the readme
p0=[t[0],0.01,0.25,1]
y_initial=exp_model(t,p0)
#plot the initial guess
plt.figure()
plt.plot(t,y,label='Data')
plt.plot(t,y_initial,label='Initial Fit')
plt.xlabel('Time')
plt.ylabel('Flux')
plt.title('Initial fit to data')
plt.legend()
plt.savefig('question2_initial_fit.pdf')

#Find best fit parameters and their error
p=newton(t,y,p0)
p_err=findError(t,y,p)
y_final=exp_model(t,p)
print('Final params are', p)
print('Errors are ', p_err)
#Plot both fits
plt.figure()
plt.plot(t,y,label='Data')
plt.plot(t,y_initial,label='Initial fit')
plt.plot(t,y_final,label='Fit minimizing Chi-square')
plt.legend()
plt.title('Initial and final fits to data')
plt.ylabel('Flux')
plt.xlabel('Time')
plt.savefig('question2_initial_and_final_fits.pdf')