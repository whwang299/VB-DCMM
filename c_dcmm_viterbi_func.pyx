# 2014 Sep16, Oct 8, Oct 16
# Wonseok Hwang
# License: GPLv3
# Tested on ipython -pylab. python 3.4 with matplotlib
# Code written by me after read
# 1. "HMM tutorial" note : http://www.ee.surrey.ac.uk/Personal/P.Jackson/tutorial/
# 2. "Sagemath hmm module" (chmm.pyx)
# 3. Double Chain Makov Model: A. Berchtold, The Double Chain Markov Model, 
#                              Technical report (Washington univ), 1999.
#. Functionalize
# 2015 Feb28, Debug: xi_arr. "denominator" hasn't been used so far.
## 0.

#from matplotlib.pylab import *
import matplotlib.pylab as plt
from os import path
from copy import deepcopy

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cdef extern from "math.h":
  double log(double)

def c_dcmm_viterbi_func( np.ndarray o_arr_filtered,
                         np.ndarray p_init_star, 
                         np.ndarray A_star, 
                         np.ndarray C_star):
                         
  """ 
  Input:
  o_arr_filtered: Noise-filtered FRET value
  Nhh: Total # of hidden-hidden states
  Nh: Total # of hiden states
  A_0: transition matrix for hh
  C_0: transition matrix for h

  Output:
  log_probability
  xhh_arr_post: Estimated one
  xh_arr: just integer-value-converted version 
  p_init_post
  A_post
  C_post
  """
  ## 1. Convert o_arr_filtered to xh_arr
  cdef int Nhh = plt.shape( A_star )[0] # or C_0[0]
  cdef int Nh = plt.shape( C_star )[1]
  cdef int Ntot = Nhh*Nh
  cdef int T = len(o_arr_filtered) # Time length

  cdef np.ndarray[DTYPE_t, ndim=1] fret_vals = np.array( list( set(o_arr_filtered) ), dtype = DTYPE)
  fret_vals.sort()

  cdef np.ndarray[DTYPE_t, ndim=1] xh_arr = np.zeros( T, dtype=DTYPE )
  cdef int i
  for t in range(T):
    #for i in range( Nh ): #len(fret_vals) ):
    for i in range( len(fret_vals) ):
      if o_arr_filtered[t] == fret_vals[i]:
        xh_arr[t] = DTYPE(i)
        continue
 
  def viterbi(xh_arr, p_init, A, C):
    """
    log-scaled viterbi
    ref) Sagemath, chmm.pyx
    """
    # init
    cdef np.ndarray[DTYPE_t, ndim=2] log_delta_arr = np.zeros( [T, Nhh], dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=2] psi = np.zeros( [T, Nhh], dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=3] log_C = np.zeros( plt.shape(C), dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=2] log_A = np.zeros( plt.shape(A), dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=1] xhh_arr_post = np.zeros( plt.shape(xh_arr), dtype=DTYPE )
    cdef int i, j, k
    cdef DTYPE_t mx 
    cdef DTYPE_t tmp
    cdef int idx 
    cdef DTYPE_t koko_max
    cdef int koko_idx

    for i in range(Nhh):
      for j in range(Nhh):
        log_A[i,j] = log( A[i,j] )

    for k in range(Nhh):
      for i in range(Nh):
        for j in range(Nh):
          log_C[k,i,j] = log( C[k,i,j] )

    for i in range(Nhh):
      log_delta_arr[1, i] = log( p_init[i] ) \
                             + log(  C[ i, int(xh_arr[0]), int(xh_arr[1]) ]  )
  
    # induction
    for t in range(2, T):
      for j in range(Nhh):
        mx = -9999999999
        idx = -1
        for i in range(Nhh):
          tmp = log_delta_arr[t-1, i] \
               + log_A[i, j]
          # to prevent -inf which causes error during cython compilation----#
          if tmp == -np.inf:                                                #
            tmp = -999999999/2.0                                            #
          # ----------------------------------------------------------------#
          if tmp > mx:
            mx = tmp
            idx = i
        log_delta_arr[t, j] = mx \
                             + log_C[ j, int(xh_arr[t-1]), int(xh_arr[t]) ]
        psi[t, j] = idx
  
    #termination
    koko_max = -9999999999
    for i in range(Nhh):
      if log_delta_arr[T-1, i] > koko_max:
        koko_max = log_delta_arr[T-1, i]
        koko_idx = i
  
    # Back tracking
    xhh_arr_post[T-1] = koko_idx
    for t in range(T-2, 0, -1): # T-2, T-1, ..., 0
      xhh_arr_post[t] = psi[t+1, int(xhh_arr_post[t+1])]
  
    #return xhh_arr_post, mx, log_delta_arr, psi
    return xhh_arr_post, log_delta_arr[-1, int(xhh_arr_post[-1])]

  return viterbi( xh_arr, p_init_star, A_star, C_star )
  
