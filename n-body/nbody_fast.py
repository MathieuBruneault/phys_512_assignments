import numpy as np
import time
import matplotlib
from bisect import bisect_left
import argparse
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm

from matplotlib  import pyplot as plt

class particles:
    def __init__(self,m=1.0,npart=1000,ngrid=100,soft=2,G=1.0,dt=0.1,size=5,BC='p',part=1):
        #Initialize some of the properties
        self.opts={}
        self.opts['soft']=soft
        self.opts['n']=npart
        self.opts['G']=G
        self.opts['dt']=dt
        self.opts['ngrid']=ngrid
        self.opts['m']=m

        #Get a random distribution of initial positions for the particles. Make sure not to start particles on the edge 
        #for non-periodic BC, since we'll be setting pot=0 there and we dont want too many particles to 'stick to the edge initially
        if BC == 'np':
            self.x=np.abs(np.random.uniform(size/ngrid,size-size/ngrid,self.opts['n']))
            self.y=np.abs(np.random.uniform(size/ngrid,size-size/ngrid,self.opts['n']))
        else:
            self.x=np.abs(np.random.uniform(0,size,self.opts['n']))
            self.y=np.abs(np.random.uniform(0,size,self.opts['n']))

        #Set some more properties (start all particles with v=0)
        self.xlim=[0, size]
        self.ylim=[0, size]
        self.vx=0*self.x
        self.vy=self.vx.copy()
        self.density_grid=np.zeros((ngrid,ngrid))
        self.BC=BC

        #In part 2 we want the particles to orbit so we start them to be close to the center, with equal but opposite velocities
        if part == 2:
            self.x[0] = 2*size/6
            self.x[1] = 4*size/6
            self.y[0] = size/2
            self.y[1] = size/2
            self.vx[0] = 0
            self.vx[1] = 0
            self.vy[0] = 1
            self.vy[1] = -1

        #The Green function will be the same for all iterations, so we compute t once
        x_grid = np.linspace(0, self.xlim[1], num=ngrid)
        y_grid = np.linspace(0, self.ylim[1], num=ngrid)
        green_fct = np.zeros((ngrid,ngrid))
        r=np.sqrt((np.array([x_grid,]*ngrid))**2+(np.array([y_grid,]*ngrid).transpose())**2)
        green_fct = 1/(4*np.pi*r)
        #Avoid infinities by setting the center of G to be G(soft) instead of G(r)
        green_fct[r < self.opts['soft']] = 1/(4*np.pi*self.opts['soft'])
        #Replicate G in all corners
        green_fct += np.flip(green_fct,1)
        green_fct += np.flip(green_fct,0)
        self.green_fct = green_fct

        if part != 4:
            self.m=np.ones(self.opts['n'])*m
        else:
            self.m=[]
            x_idx = [take_closest(x_grid,i) for i in self.x]
            y_idx = [take_closest(y_grid,i) for i in self.y]
            r=np.sqrt((np.array([x_grid,]*ngrid)-size/2)**2+(np.array([y_grid,]*ngrid).transpose()-size/2)**2)
            for(i,j) in zip(x_idx,y_idx):
                self.m=np.append(self.m,m/r[i][j]**3)
            self.m[self.m > m/(self.opts['soft'])**3] = m/(self.opts['soft'])**3

    def get_density(self):
        #Create a grid to discretize position
        grid = np.zeros((self.opts['ngrid'],self.opts['ngrid']))
        x_grid = np.linspace(self.xlim[0], self.xlim[1], num=self.opts['ngrid'])
        y_grid = np.linspace(self.ylim[0], self.ylim[1], num=self.opts['ngrid'])
        #Find in which grid cell each particle belongs
        x_idx = [take_closest(x_grid,i) for i in self.x]
        y_idx = [take_closest(y_grid,i) for i in self.y]
        k=0
        #Loop over all particles and associate a mass to each grid cell in which there are particles
        #and divide by the grid spacing squared to get density
        for (i,j) in zip(x_idx,y_idx):
            grid[i][j] += self.m[k]
            k+=1
        grid /= ((x_grid[1]-x_grid[0])*(y_grid[1]-y_grid[0]))
        self.density_grid=grid
        return x_idx, y_idx

    def get_forces(self):
        self.fx=np.zeros(self.opts['n'])
        self.fy=0*self.fx
        #Update the density grid and get the corresponding of particles in that grid
        x_idx, y_idx = part.get_density()
        #Convolve density grid and Gto get the potential
        pot = np.fft.ifft2(np.fft.fft2(self.density_grid)*np.fft.fft2(self.green_fct)).real
        pot=0.5*(np.roll(pot, 1, axis=1) + pot)
        pot=0.5*(np.roll(pot, 1, axis=0) + pot)
        #Set the potential to 0 on the edges for on periodic BC
        if BC == 'np':
            pot[0,:]=0
            pot[-1,:]=0
            pot[:,0]=0
            pot[:,-1]=0
        #Numerically evaluate forces on the grid and associate them to their corresponding particles
        fy = -0.5*(np.roll(pot, 1, axis = 1) - np.roll(pot, -1, axis=1))*self.density_grid*(self.xlim[1]/self.opts['ngrid'])
        fx = -0.5*(np.roll(pot, 1, axis = 0) - np.roll(pot, -1, axis=0))*self.density_grid*(self.ylim[1]/self.opts['ngrid'])
        for (i,j) in enumerate(x_idx):
            self.fx[i] = fx[j][y_idx[i]]
            self.fy[i] = fy[j][y_idx[i]]
        return -0.5*np.sum(pot)

    def evolve(self):
        #Evolve the positions base on previously calculated velocities
        self.x+=self.vx*self.opts['dt']
        self.y+=self.vy*self.opts['dt']
        #We want particles to be on a torus for periodic BC, so they loop if they go too far
        if self.BC == 'p':
            self.x = self.x%(self.xlim[1])
            self.y = self.y%(self.ylim[1])
        #Get potential energy and update forces. Update velocities based on forces
        pot=self.get_forces()
        self.vx+=self.fx*self.opts['dt']
        self.vy+=self.fy*self.opts['dt']
        #Compute kinetic energy and return total energy
        kinetic=0.5*np.sum(self.m*(self.vx**2+self.vy**2))
        return pot+kinetic


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns index of closest value to myNumber.

    If two numbers are equally close, return index of the smallest number.

    Makes an O(n) process O(log(n))
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return len(myList)-1
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1

if __name__=='__main__':
    #Get arguments to know which part to run and with which boundary conditions
    arg_parser = argparse.ArgumentParser(
        description="Numerically solve the n-body problem"
    )
    arg_parser.add_argument("part_number", help="Which part should be executed (1, 2, 3 or 4 are accepted)", type=int)
    arg_parser.add_argument("-BC", help="Boundary conditions to use. Accepted values are p (periodic) and np (non periodic). Assumes periodic BC if none is given", type=str, default='p')
    args = arg_parser.parse_args()

    #Check that the BC argument is valid
    if args.BC == 'p' or args.BC == 'np':
        BC = args.BC
    else:
        print('Not a valid boundary condition, use -h for more information')
        exit()

    #Initialize the problem differently depending on which part we re running (and check that the part number is valid)
    if args.part_number == 1:
        n=1
        ngrid=50
        m=1.0/n
        dt=1
        soft=1
    elif args.part_number == 2:
        n=2
        ngrid=100
        m=1.0/n
        dt=0.03
        soft=1
    elif args.part_number == 3:
        n=100000
        ngrid=512
        m=1/n
        dt=0.1
        soft=0.1
    elif args.part_number == 4:
        n=100000
        ngrid=512
        m=1/n
        dt=0.001
        soft=0.1
    else:
        print('Not a valid part number, use -h for more information')
        exit()

    #Intialize the particles class and get a particle object
    part=particles(m=m,npart=n,ngrid=ngrid,dt=dt,soft=soft,BC=BC,part=args.part_number)
    grid=[]
    x=[]
    y=[]
    all_energy=[]

    #Loop through time and evolve the particles each iteration
    for i in range(0,750):
        energy=part.evolve()
        all_energy.append(energy)
        print('Energy is ', energy)
        plt.clf()
        if args.part_number == 1 or args.part_number == 2:
            #Plot the particles directly if part 1 or 2
            plt.scatter(part.x,part.y)
            plt.ylim([0,5])
            plt.xlim([0,5])
            plt.title('Part ' + str(args.part_number) + ' with ' + str(BC) + 'eriodic boundary conditions')
            plt.xlabel('x')
            plt.ylabel('y')
            x.append(part.x)
            y.append(part.y)
        else:
            #Plot the particle density if we are doing part 3 or 4
            plt.pcolormesh(part.density_grid, norm=LogNorm())
            plt.colorbar()
            plt.title('Part ' + str(args.part_number) + ' with ' + str(BC) + 'eriodic boundary conditions')
            plt.xlabel('x')
            plt.ylabel('y')
            grid.append(part.density_grid)
        plt.pause(1e-3)
    if args.part_number == 1 or args.part_number == 2:
        fig,ax = plt.subplots(figsize=(8,6))
        cax=ax.scatter(x[0],y[0])
        ax.set_title('Part ' + str(args.part_number) + ' with ' + str(BC) + 'eriodic boundary conditions')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim([0,5])
        ax.set_ylim([0,5])
        def animate(i):
            cax=ax.scatter(x[i],y[i])
        anim=FuncAnimation(fig,animate,interval=40,frames=750,repeat=True,blit=False,save_count=1000)
        anim.save('Part_'+str(args.part_number)+'_'+BC+'.mp4')
    else:
        fig,ax = plt.subplots(figsize=(8,6))
        cax=ax.imshow(grid[0])
        ax.set_title('Part ' + str(args.part_number) + ' with ' + str(BC) + 'eriodic boundary conditions')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        def animate(i):
            cax.set_array(grid[i])
        anim=FuncAnimation(fig,animate,interval=40,frames=750,repeat=True,blit=False,save_count=1000)
        anim.save('Part_'+str(args.part_number)+'_'+BC+'.mp4')

    plt.figure()
    plt.plot(np.arange(len(all_energy)),all_energy)
    plt.xlabel('Iteration')
    plt.ylabel('Energy')
    plt.title('Part ' + str(args.part_number) + ' with ' + str(BC) + 'eriodic boundary conditions - Energy')
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.savefig('Part_'+str(args.part_number)+'_'+BC+'_energy.pdf')