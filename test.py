
# The random number generator is the same in both the C
# and python version, since its taken from gsl.

from ctypes import *
import pygsl
import pygsl.rng
import numpy as np
from smm import *

r = pygsl.rng.ranlxs1()
seed = 1214314
r.set(seed)

c = 0.12
k = (2 * np.pi)/7.5
T = (2 * np.pi)/(c*k)
randomshift = r.uniform() * T;

print "Seed = ", seed, "RandomShift = ", randomshift

Nsensors = 200
StreamFunction = 5 # Meandering Jet
Params = (c_double*10)(1.2, c, k, 0.4, 0.3, randomshift )
timestep = 0.01
SimStep = timestep*0.03*24*60*60
SMM_init(Nsensors, StreamFunction, Params, 1, timestep)
NDeploy = Nsensors
SMM_deploy_nodes(NDeploy)

x = np.zeros(Nsensors, dtype=np.float64)
y = np.zeros(Nsensors, dtype=np.float64)

for i in range(Nsensors):
	x[i] = r.uniform() * 4
	y[i] = r.uniform() * 4 - 2

print "%.2f\t%.3f\t%.3f" % (0,min(x),max(x))

for i in range(100):
	MoveSensors(x,y)
	print "%.2f\t%.2f\t%.2f" % ((i+1)*SimStep,min(x)*1000,max(x)*1000)



