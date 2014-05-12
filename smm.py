import numpy as np
from ctypes import *


# load the library

smm = cdll.LoadLibrary("./libsmm.so")
smm.SMM_init.argtypes = [c_int, c_int, POINTER(c_double), c_int, c_double]
smm.SMM_init.restype = None
smm.rk2.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C,W'),
                    np.ctypeslib.ndpointer(dtype=np.float64,ndim=1,flags='C,W')]
smm.rk2.restype = None

MoveSensors = smm.rk2
SMM_init = smm.SMM_init
SMM_deploy_nodes = smm.SMM_deploy_nodes
