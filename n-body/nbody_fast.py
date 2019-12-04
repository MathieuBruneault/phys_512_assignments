import numpy as np
import time
import matplotlib
from bisect import bisect_left
from scipy.signal import convolve2d
import argparse

from matplotlib  import pyplot as plt

class particles:
    def __init__(self,m=1.0,npart=1000,ngrid=100,soft=2,G=1.0,dt=0.1,size=5,BC='p',part=1):
        self.opts={}
        self.opts['soft']=soft
        self.opts['n']=npart
        self.opts['G']=G
        self.opts['dt']=dt
        self.opts['ngrid']=ngrid
        self.opts['m']=m

        self.x=np.abs(np.random.uniform(0,size,self.opts['n']))
        self.y=np.abs(np.random.uniform(0,size,self.opts['n']))
        self.xlim=[0, size]
        self.ylim=[0, size]
        self.m=np.ones(self.opts['n'])*m
        self.vx=0*self.x
        self.vy=self.vx.copy()
        self.density_grid=np.zeros((ngrid,ngrid))
        self.BC=BC

        if part == 2:
            self.x[0] = 2*size/6
            self.x[1] = 4*size/6
            self.y[0] = size/2
            self.y[1] = size/2
            self.vx[0] = 0
            self.vx[1] = 0
            self.vy[0] = 1
            self.vy[1] = -1

        x_grid = np.linspace(0, self.xlim[1], num=ngrid)
        y_grid = np.linspace(0, self.ylim[1], num=ngrid)
        green_fct = np.zeros((ngrid,ngrid))
        #r=np.sqrt((np.array([x_grid,]*ngrid)-x_grid[ngrid//2])**2+(np.array([y_grid,]*ngrid).transpose()-y_grid[ngrid//2])**2)
        r=np.sqrt((np.array([x_grid,]*ngrid))**2+(np.array([y_grid,]*ngrid).transpose())**2)
        green_fct = 1/(4*np.pi*r)
        green_fct[r < self.opts['soft']] = 1/(4*np.pi*self.opts['soft'])
        #green_fct += 1/(4*np.pi*self.opts['soft'])
        green_fct += np.flip(green_fct,1)
        green_fct += np.flip(green_fct,0)
        self.green_fct = green_fct

    def get_density(self):
        grid = np.zeros((self.opts['ngrid'],self.opts['ngrid']))
        x_grid = np.linspace(self.xlim[0], self.xlim[1], num=self.opts['ngrid'])
        y_grid = np.linspace(self.ylim[0], self.ylim[1], num=self.opts['ngrid'])
        x_idx = [take_closest(x_grid,i) for i in self.x]
        y_idx = [take_closest(y_grid,i) for i in self.y]
        for (i,j) in zip(x_idx,y_idx):
            grid[i][j] += self.opts['m']
        grid /= ((x_grid[1]-x_grid[0])*(y_grid[1]-y_grid[0]))
        self.density_grid=grid
        return x_idx, y_idx

    def get_forces(self):
        self.fx=np.zeros(self.opts['n'])
        self.fy=0*self.fx
        x_idx, y_idx = part.get_density()
        pot = np.fft.ifft2(np.fft.fft2(self.density_grid)*np.fft.fft2(self.green_fct)).real
        pot=0.5*(np.roll(pot, 1, axis=1) + pot)
        pot=0.5*(np.roll(pot, 1, axis=0) + pot)
        if self.BC == 'p':
            fy = -0.5*(np.roll(pot, 1, axis = 1) - np.roll(pot, -1, axis=1))*self.density_grid*(self.xlim[1]/self.opts['ngrid'])
            fx = -0.5*(np.roll(pot, 1, axis = 0) - np.roll(pot, -1, axis=0))*self.density_grid*(self.ylim[1]/self.opts['ngrid'])
        for (i,j) in enumerate(x_idx):
            self.fx[i] = fx[j][y_idx[i]]
            self.fy[i] = fy[j][y_idx[i]]
        return -0.5*np.sum(pot)

    def evolve(self):
        self.x+=self.vx*self.opts['dt']
        self.y+=self.vy*self.opts['dt']
        if self.BC == 'p':
            self.x = self.x%(self.xlim[1])
            self.y = self.y%(self.ylim[1])
        pot=self.get_forces()
        self.vx+=self.fx*self.opts['dt']
        self.vy+=self.fy*self.opts['dt']
        kinetic=0.5*np.sum(self.m*(self.vx**2+self.vy**2))
        return pot+kinetic


def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return len(myList)
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return pos
    else:
       return pos - 1

if __name__=='__main__':
    arg_parser = argparse.ArgumentParser(
        description="Numerically solve the n-body problem"
    )
    arg_parser.add_argument("part_number", help="Which part should be executed (1, 2, 3 or 4 are accepted)", type=int)
    arg_parser.add_argument("-BC", help="Boundary conditions to use. Accepted values are p (periodic) and np (non periodic). Assumes periodic BC if none is given", type=str, default='p')
    args = arg_parser.parse_args()
    if args.BC == 'p' or args.BC == 'np':
        BC = args.BC
    else:
        print('Not a valid boundary condition, use -h for more information')
        exit()
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
        dt=5
        soft=2
    else:
        print('Not a valid part number, use -h for more information')
        exit()

    part=particles(m=m,npart=n,ngrid=ngrid,dt=dt,soft=soft,BC=BC,part=args.part_number)

    for i in range(0,10000):
        energy=part.evolve()
        print('Energy is ', energy)
        plt.clf()
        if args.part_number == 1 or args.part_number == 2:
            plt.scatter(part.x,part.y)
            plt.ylim([0,5])
            plt.xlim([0,5])
            plt.title('Part ' + str(args.part_number) + ' with ' + str(BC) + 'eriodic boundary conditions')
            plt.xlabel('x')
            plt.ylabel('y')
        else:
            plt.pcolormesh(part.density_grid)
            plt.colorbar()
            plt.title('Part ' + str(args.part_number) + ' with ' + str(BC) + 'eriodic boundary conditions')
            plt.xlabel('x')
            plt.ylabel('y')
        plt.pause(1e-3)