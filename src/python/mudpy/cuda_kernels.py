'''
Marc Coll 05/2021

CUDA kernels for those methods that use the GPU
'''

import numpy as np
import math
import numba
from numba import cuda


@cuda.jit
def tshift_trace_kernel(nss,ess,zss,nds,eds,zds,
                        Nss_data_array,Ess_data_array,Zss_data_array,Nds_data_array,Eds_data_array,Zds_data_array,
                        index,npts,nshift):
    '''
    Does pretty much the same as tshift_trace, but more efficiently and adapted to run on the GPU
    Warning: nshift has to be calculated outside of the kernel. Otherwise, it crashes. Probably
    due to some weird compilation bug
    '''

    # Identify the thread in order to know which one of the array elements is ours
    pos = cuda.grid(1)

    if (pos < npts):
        if (nshift > 0):
            #This code does the equivalent to:
            #   arr.data=r_[zeros(nshift),arr.data[0:npts-nshift]]
            if (pos < nshift):
                nss[pos]=0.0
                ess[pos]=0.0
                zss[pos]=0.0
                nds[pos]=0.0
                eds[pos]=0.0
                zds[pos]=0.0
            else:
                nss[pos]=Nss_data_array[index][pos-nshift]
                ess[pos]=Ess_data_array[index][pos-nshift]
                zss[pos]=Zss_data_array[index][pos-nshift]
                nds[pos]=Nds_data_array[index][pos-nshift]
                eds[pos]=Eds_data_array[index][pos-nshift]
                zds[pos]=Zds_data_array[index][pos-nshift]
        else:
            #This code does the equivalent to:
            #   arr.data=r_[arr.data[-nshift:],arr.data[-1]*ones(-nshift)]
            if (pos < (npts + nshift)):
                nss[pos]=Nss_data_array[index][pos-nshift]
                ess[pos]=Ess_data_array[index][pos-nshift]
                zss[pos]=Zss_data_array[index][pos-nshift]
                nds[pos]=Nds_data_array[index][pos-nshift]
                eds[pos]=Eds_data_array[index][pos-nshift]
                zds[pos]=Zds_data_array[index][pos-nshift]
            else:
                nss[pos]=Nss_data_array[index][-1]
                ess[pos]=Ess_data_array[index][-1]
                zss[pos]=Zss_data_array[index][-1]
                nds[pos]=Nds_data_array[index][-1]
                eds[pos]=Eds_data_array[index][-1]
                zds[pos]=Zds_data_array[index][-1]


@cuda.jit
def triangle_stf_kernel(Mdot,t,rise_time,dt,time_offset,m,b1,b2):
    # Identify the thread in order to know which one of the array elements is ours
    i = cuda.grid(1)

    if (i < len(Mdot)):
        t[i] = dt * i

        if (t[i] <= rise_time / 2):
            Mdot[i] = m * t[i] + b1
        elif (t[i] > rise_time):
            Mdot[i] = 0.0
        elif (t[i] > rise_time / 2):
            Mdot[i] = -m * t[i] + b2
        else:
            Mdot[i] = 0.0

        #offset origin time
        t[i] = t[i] + time_offset


@cuda.jit
def cosine_stf_kernel(Mdot,t,rise_time,dt,time_offset,tau1,tau2,Cn):
    # Identify the thread in order to know which one of the array elements is ours
    i = cuda.grid(1)

    if (i < len(Mdot)):
        t[i] = dt * i

        if (t[i] < tau1):
            Mdot[i] = Cn * (0.7 - 0.7 * math.cos(t[i] * math.pi / tau1) +0.6 * math.sin(0.5 * math.pi * t[i] / tau1))
        elif ((t[i] >= tau1) & (t[i] < 2 * tau1)):
            Mdot[i] = Cn * (1.0 - 0.7 * math.cos(t[i] * math.pi / tau1) + 0.3 * math.cos(math.pi * (t[i] - tau1) / tau2))
        elif ((t[i] >= 2 * tau1) & (t[i] < rise_time)):
            Mdot[i] = Cn * (0.3 + 0.3 * math.cos(math.pi * (t[i] - tau1) / tau2))
        else:
            Mdot[i] = 0.0

        #offset origin time
        t[i] = t[i] + time_offset


@cuda.jit
def dreger_stf_kernel(Mdot,t,dt,time_offset,zeta,tau):
    # Identify the thread in order to know which one of the array elements is ours
    i = cuda.grid(1)

    if (i < len(Mdot)):
        t[i] = dt * i

        Mdot[i] = (t[i] ** zeta) * math.exp(-t[i] / tau)

        #offset origin time
        t[i] = t[i] + time_offset
