#!/usr/bin/env/python
import time
from mpi4py import MPI
import numpy as np
Ntot = 1e7

comm=MPI.COMM_WORLD
myrank = comm.Get_rank()
nproc=comm.Get_size()
tstart = time.clock()
N=int(Ntot)/nproc
pi = 0
for i in range(N): #N is the local N
	x = np.random.rand()
	y = np.random.rand()
	r = np.sqrt(x**2+y**2)
	if r<1:
		pi+=1
pi = comm.reduce(pi,op=MPI.SUM,root=0)
Ntot = comm.reduce(N,op=MPI.SUM,root=0)
tstop = time.clock()
if myrank==0:
	pi*=4.0/Ntot
	relerr = abs(pi-np.pi)/np.pi
print "approx pi : {}".format(pi)
print "relative error :{}".format(relerr)
print "Elapsed time :{}".format(tstop-tstart)