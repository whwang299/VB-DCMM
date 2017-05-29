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
# 2015 Feb27, Variational Bayes version
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

def c_dcmm4vb_func( np.ndarray o_arr_filtered,
                 np.ndarray p_init_0, 
                 np.ndarray A_0, 
                 np.ndarray C_0):
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
  cdef int Nhh = plt.shape( A_0 )[0] # or C_0[0]
  cdef int Nh = plt.shape( C_0 )[1]
  cdef int Ntot = Nhh*Nh
  cdef int T = len(o_arr_filtered) # Time length of observable
  cdef int Th = len(o_arr_filtered) - 1 # Time length of hidden state

  cdef np.ndarray[DTYPE_t, ndim=1] fret_vals = np.array( list( set(o_arr_filtered) ), dtype = DTYPE)
  fret_vals.sort()
  cdef np.ndarray[DTYPE_t, ndim=1] xh_arr = np.zeros( T, dtype=DTYPE )
  cdef int i
  cdef DTYPE_t val
  for t in range(T):
    #for i in range( Nh ): #len(fret_vals) ):
    for i in range( len(fret_vals) ):
      if o_arr_filtered[t] == fret_vals[i]:
        xh_arr[t] = DTYPE(i)
        continue
  
  def construct_scaled_alpha(np.ndarray[DTYPE_t, ndim=1] xh_arr,
                             np.ndarray[DTYPE_t, ndim=1] p_init, 
                             np.ndarray[DTYPE_t, ndim=2] A, 
                             np.ndarray[DTYPE_t, ndim=3] C):
    """
    up to T-1
    xh(1) determined by xh(0) and xhh(1).
    alpha_arr[0]: blank. Just 0.
    a_scale[0]: blank
    """
    cdef np.ndarray[DTYPE_t, ndim=2] alpha_arr_r = np.zeros( [T, Nhh], dtype=DTYPE) # r stands for rescaled
    cdef np.ndarray [DTYPE_t, ndim=1] a_scale = np.zeros( T, dtype=DTYPE )
    #zeros([T,1]) --> 2-dim array, zeros(T) --> 1-dim array
    cdef DTYPE_t log_probability
  
    cdef int i, j, k, t
    for i in range(Nhh):
      alpha_arr_r[1,i] = p_init[i] * C[i, int(xh_arr[0]), int(xh_arr[1]) ]
      a_scale[1] += alpha_arr_r[1,i]
  
    for i in range(Nhh):
      alpha_arr_r[1,i] /= a_scale[1]
  
    for t in range(2, T):
      for j in range(Nhh):
        for i in range(Nhh):
          alpha_arr_r[t, j] += alpha_arr_r[t-1, i] * A[i, j]

        alpha_arr_r[t,j] *= C[ j, int(xh_arr[t-1]), int(xh_arr[t]) ]
        a_scale[t] += alpha_arr_r[t, j]
  
      for j in range(Nhh):
        alpha_arr_r[t,j] /= a_scale[t]
  

    for t in range(1,T):
      log_probability += log( a_scale[t] )
    return alpha_arr_r, a_scale, log_probability
  
  
  def construct_scaled_beta(np.ndarray[DTYPE_t, ndim=1] xh_arr,
                            np.ndarray[DTYPE_t, ndim=1] p_init, 
                            np.ndarray[DTYPE_t, ndim=2] A, 
                            np.ndarray[DTYPE_t, ndim=3] C, 
                            np.ndarray[DTYPE_t, ndim=1] a_scale):
    """
    beta_arr[0]: blank
    """
    cdef np.ndarray[DTYPE_t, ndim=2] beta_arr_r = np.zeros( [T, Nhh], dtype=DTYPE )
    cdef int i, j, t

    for i in range(Nhh):
      beta_arr_r[ T-1, i] = 1 # arr[T-1] indicates T th value

    for t in range(T-2, 0, -1): # T-2, T-3, ..., 1
      for i in range(Nhh):
        for j in range(Nhh):
          beta_arr_r[t,i] += C[ j, int(xh_arr[t]), int(xh_arr[t+1]) ] \
                            * A[i, j] \
                            * beta_arr_r[t+1, j]
        
      for i in range(Nhh):
        beta_arr_r[t, i] /= a_scale[t]
  
    return beta_arr_r
  
  def construct_gamma(np.ndarray[DTYPE_t, ndim=2] alpha_arr_r, 
                      np.ndarray[DTYPE_t, ndim=2]  beta_arr_r):
    """
    gamma_arr[0]: blank
    """
    cdef int t
    cdef int j
    cdef int k
    cdef np.ndarray[DTYPE_t, ndim=2] gamma_arr = np.zeros( [T, Nhh], dtype=DTYPE)
    cdef DTYPE_t denominator

    for t in range(1, T):
      denominator = 0.0
      for j in range(Nhh):
        gamma_arr[t, j] = alpha_arr_r[t, j] * beta_arr_r[t, j]
        denominator += gamma_arr[t, j]
  
      if denominator > 0:
        for k in range(Nhh):
          gamma_arr[t,k] /= denominator
      else: # prevent nan value of gamma 
        for k in range(Nhh):
          gamma_arr[t,k] = 1.0/Nhh
  
    return gamma_arr
  
  def construct_xi( np.ndarray[DTYPE_t, ndim=2] alpha_arr_r, 
                    np.ndarray[DTYPE_t, ndim=2] beta_arr_r, 
                    np.ndarray[DTYPE_t, ndim=1] xh_arr,
                    np.ndarray[DTYPE_t, ndim=2] A, 
                    np.ndarray[DTYPE_t, ndim=3] C ):
    """
    xi_arr[t=0]: blank
    xi_arr[t=T-1]: also blank (as xi_arr requires t and t+1
      """
    cdef int t, i, j
    cdef np.ndarray[DTYPE_t, ndim=3] xi_arr = np.zeros( [T, Nhh, Nhh], dtype=DTYPE)
    cdef DTYPE_t denominator

    for t in range(1, T-1): 
      denominator = 0.0
      for i in range(Nhh):
        for j in range(Nhh):
          xi_arr[t, i, j] = alpha_arr_r[t, i] \
                         * A[i, j] \
                         * C[j, int(xh_arr[t]), int(xh_arr[t+1])] \
                         * beta_arr_r[t+1, j]
          denominator += xi_arr[t, i, j]
  
      if denominator < 0:
        for i in range(Nhh):
          for j in range(Nhh):
            xi_arr[t, i, j] = 1.0 / Nhh / Nhh
  
    return xi_arr/denominator
  
  def get_gamma_xi_logL( np.ndarray[DTYPE_t, ndim=1] xh_arr,
                    np.ndarray[DTYPE_t, ndim=1] p_init_in, 
                    np.ndarray[DTYPE_t, ndim=2] A_in, 
                    np.ndarray[DTYPE_t, ndim=3] C_in): 

    alpha_arr_r, a_scale, log_probability = construct_scaled_alpha(
                                              xh_arr,
                                              p_init_in,
                                              A_in,
                                              C_in)
    beta_arr_r = construct_scaled_beta( 
                    xh_arr,
                    p_init_in,
                    A_in,
                    C_in,
                    a_scale)

    gamma_arr = construct_gamma(alpha_arr_r, beta_arr_r)

    xi_arr = construct_xi( alpha_arr_r, 
                           beta_arr_r,
                           xh_arr,
                           A_in,
                           C_in )
  
 
    return gamma_arr, xi_arr, log_probability

  return get_gamma_xi_logL( o_arr_filtered, p_init_0, A_0, C_0 )
                              
