import matplotlib.pyplot as plt
import numpy as np

chains=np.loadtxt('q4_chains.txt')
nsteps=len(chains[:,0])
x=np.arange(nsteps)

print('Params are ', np.mean(chains, axis=0))
cov=np.cov(chains[int(nsteps/2):,:].T)
err=np.sqrt(np.diag(cov))
print('Errors are ', err)

fig, axs = plt.subplots(6)
axs[0].plot(x, chains[:,0])
axs[0].set(ylabel='H0')
axs[1].plot(x, chains[:,1])
axs[1].set(ylabel='wb')
axs[2].plot(x, chains[:,2])
axs[2].set(ylabel='wc')
axs[3].plot(x, chains[:,3])
axs[3].set(ylabel='tau')
axs[4].plot(x, chains[:,4])
axs[4].set(ylabel='As')
axs[5].plot(x, chains[:,5])
axs[5].set(ylabel='slope')
axs[5].set(xlabel='Iteration')
plt.show()

