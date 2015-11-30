# 2014 Sep16, Oct 8
# Wonseok Hwang
# License: GPLv3
# Tested on ipython -pylab. python 3.4 with matplotlib
# Code is written by Wonseok Hwang after reading
# 1. "HMM tutorial" note : http://www.ee.surrey.ac.uk/Personal/P.Jackson/tutorial/
# 2. "Sagemath hmm module" (chmm.pyx)
## 0.
from matplotlib.pylab import *
from os import path

import numpy as np
cimport numpy as np

DTYPE = np.float64

ctypedef np.float64_t DTYPE_t
# Now following functions are overloaded by C-function
cdef extern from "math.h":
  double sqrt(double)
  double exp(double)
  double log(double)
  double fabs(double)

## 1. Parameters
sig_min = 0.001

def c_hmm_func(o_arr, p_init_0, p_xh_0, mu_arr_0, sig_arr_0):
  """
  input:
  o_arr
  p_init_0
  p_xh_0
  mu_arr_0
  sig_arr_0

  output:
  log_probability
  xh_arr_post
  p_xh_post
  mu_arr_post
  sig_arr_post
  """
  ## 2. Load data (generated from Matlab/Octave file)
  cdef int T = len(o_arr)
  cdef int N = len(p_xh_0)

  def emmission_prob(o, i, mu, sig):
    """
    o: Observable
    i: hidden state
    mu_arr: means
    sig_arr: stds
    """
    cdef DTYPE_t c_pi = pi
    cdef DTYPE_t c_o = o
    cdef DTYPE_t c_mu = mu
    cdef DTYPE_t c_sig = sig

    return 1.0 / sqrt( 2 * c_pi ) / c_sig \
        * exp( -(c_o - c_mu)**2 / ( 2 * c_sig**2 ) )
  
  def construct_scaled_alpha(np.ndarray[DTYPE_t, ndim=1] o_arr, 
                             np.ndarray[DTYPE_t, ndim=1] p_init,
                             np.ndarray[DTYPE_t, ndim=2] p_xh,
                             np.ndarray[DTYPE_t, ndim=1] mu_arr, 
                             np.ndarray[DTYPE_t, ndim=1] sig_arr):
    """
    Forward algorithm with scaling
    construct alpha array on give o_arr and parameters
    a(t,i) = P(o1, o2, ...,ot, xt = i | lambda)
    o_arr: observation array
    p_init: pi (initial probability)
    p_xh: transition matrix
    mu_o: means of emmision of each state
    sig_o: stds of emmision of each state
  
  
    alpha_r[t4] = alpha_ori[t4] / (s4 * s3 * s2 * s1)
    where s4 = sum( 1/s3/s3/s1 * alpha_arr_ori[t4] )
    s3 = sum( alpha_r[t3, i], i )
    a_scale[t4] = s4 
  
    log(s1) + log(s2) + log(s3) + ... + log(sT)
     =  log (s1 * ... * sT-1  * 1/(s1* ... * sT-1) * sum( alpha_r[T, i], i)
     = log ( P(obs | lambda) )
  
    return
    a_scale: alpha scaler NOT w.r.t. real alpha but will be used for beta calculation
  
    """
    cdef np.ndarray[DTYPE_t, ndim=2] alpha_arr_r = np.zeros( [T, N], dtype=DTYPE ) # r stands for rescaled
    cdef np.ndarray[DTYPE_t, ndim=1] a_scale = zeros( T, dtype=DTYPE )
    cdef int i
    cdef int j
    cdef int t
    cdef DTYPE_t c_log_probability = 0.0
  
    for i in range(N):
      alpha_arr_r[0,i] = p_init[i] * emmission_prob(o_arr[0], i, mu_arr[i], sig_arr[i])
      a_scale[0] += alpha_arr_r[0,i]
  
    #alpha_arr_r[0,:] = alpha_arr_r[0,:] / a_scale[0]
    for j in range(N):
      alpha_arr_r[0,j] /= a_scale[0] # /= is much faster then = ## / ## notation
  
    for t in range(1, T):
      for j in range(N):
        for i in range(N):
          alpha_arr_r[t, j] +=  alpha_arr_r[t-1, i] * p_xh[i, j] 

        alpha_arr_r[t,j] *= emmission_prob( o_arr[t], j, mu_arr[j], sig_arr[j] )
        a_scale[t] += alpha_arr_r[t, j]

      for j in range(N):
        alpha_arr_r[t,j] /= a_scale[t]
  
    for t in range(T):
      c_log_probability += log(a_scale[t])
    return alpha_arr_r, a_scale, c_log_probability
    #return alpha_arr_r, a_scale, sum( log( a_scale ) )
     
  
  def construct_scaled_beta(np.ndarray[DTYPE_t, ndim=1] o_arr, 
                            np.ndarray[DTYPE_t, ndim=1] p_init,
                            np.ndarray[DTYPE_t, ndim=2] p_xh, 
                            np.ndarray[DTYPE_t, ndim=1] mu_arr, 
                            np.ndarray[DTYPE_t, ndim=1] sig_arr, 
                            np.ndarray[DTYPE_t, ndim=1] a_scale):
    """
    Forward algorithm with scaling
    construct alpha array on give o_arr and parameters
    a(t,i) = P(o1, o2, ...,ot, xt = i | lambda)
    o_arr: observation array
    p_init: pi (initial probability)
    p_xh: transition matrix
    mu_o: means of emmision of each state
    sig_o: stds of emmision of each state
  
  
    alpha_r[t4] = alpha_ori[t4] / (s4 * s3 * s2 * s1)
    where s4 = sum( 1/s3/s3/s1 * alpha_arr_ori[t4] )
    s3 = sum( alpha_r[t3, i], i )
    a_scale[t4] = s4 
  
    log(s1) + log(s2) + log(s3) + ... + log(sT)
     =  log (s1 * ... * sT-1  * 1/(s1* ... * sT-1) * sum( alpha_r[T, i], i)
     = log ( P(obs | lambda) )
  
    return
    a_scale: alpha scaler NOT w.r.t. real alpha but will be used for beta calculation
  
    """
    cdef np.ndarray[DTYPE_t, ndim=2] beta_arr_r = np.zeros( [T, N], dtype=DTYPE )
    cdef DTYPE_t c_b_sum=0
    cdef int j
    cdef int t
    cdef int i

    for j in range(N):
      beta_arr_r[ T-1, j] = 1 / a_scale[T-1] # arr[T-1] indicates T th value
    for t in range(T-2, -1, -1): # T-2, T-1, ..., 0
      for i in range(N):
        c_b_sum = 0
        for j in range(N):
          c_b_sum += p_xh[i, j] \
                    * emmission_prob( o_arr[t+1], j, mu_arr[j], sig_arr[j] ) \
                    * beta_arr_r[t+1, j]
  
        beta_arr_r[t, i] = c_b_sum / a_scale[t]
    
    return beta_arr_r
  
  def construct_gamma(np.ndarray[DTYPE_t, ndim=2] alpha_arr_r, 
                      np.ndarray[DTYPE_t, ndim=2] beta_arr_r):
    cdef np.ndarray[DTYPE_t, ndim=2] gamma_arr = np.zeros( [T, N], dtype=DTYPE )
    cdef DTYPE_t c_denom
    cdef int t
    cdef int j
    cdef int k
    for t in range(T):
  
      c_denom = 0.0
      #denominator = 0.0
      for j in range(N):
        gamma_arr[t, j] = alpha_arr_r[t, j] * beta_arr_r[t, j]
        c_denom += gamma_arr[t, j]
        #denominator += gamma_arr[t, j]
  
      if c_denom > 0:
        for k in range(N):
          gamma_arr[t,k] /= c_denom
      else: # prevent nan value of gamma 
        for k in range(N):
          gamma_arr[t,k] = 1.0/N
  
    return gamma_arr
  
  def construct_xi( np.ndarray[DTYPE_t, ndim=2] alpha_arr_r, 
                    np.ndarray[DTYPE_t, ndim=2] beta_arr_r, 
                    np.ndarray[DTYPE_t, ndim=1] o_arr, 
                    np.ndarray[DTYPE_t, ndim=2] p_xh, 
                    np.ndarray[DTYPE_t, ndim=1] mu_arr, 
                    np.ndarray[DTYPE_t, ndim=1] sig_arr ):
    """
    For xi, Sage convention used: t=[1, T-1]
    xi(t, i, j) ~ alpha(t) * Aij * bj( o(t+1) ) * beta(t+1)
    
    (different with ppt convention: t= [2, T])
    xi_ppt(t, i, j) ~ alpha(t-1) * Aij * bj(ot) * beta(t)
    """

    cdef np.ndarray[DTYPE_t, ndim=3] xi_arr = np.zeros( [T, N, N], dtype=DTYPE )
    cdef DTYPE_t c_denom
    cdef int t
    cdef int i
    cdef int j
    for t in range(T-1):
      c_denom = 0.0
      #denominator = 0.0
      for i in range(N):
        for j in range(N):
          xi_arr[t, i, j] = alpha_arr_r[t, i] \
                         * p_xh[i, j] \
                         * emmission_prob( o_arr[t+1], j, mu_arr[j], sig_arr[j]) \
                         * beta_arr_r[t+1, j]
          c_denom += xi_arr[t, i, j]
  
      if c_denom < 0:
        for i in range(N):
          for j in range(N):
            xi_arr[t, i, j] = 1.0 / N / N
  
    return xi_arr
  
  def construct_params(np.ndarray[DTYPE_t, ndim=2] gamma_arr, 
                       np.ndarray[DTYPE_t, ndim=3] xi_arr, 
                       np.ndarray[DTYPE_t, ndim=1] o_arr, 
                       np.ndarray[DTYPE_t, ndim=1] mu_arr_in, 
                       np.ndarray[DTYPE_t, ndim=1] sig_arr_in):
    cdef np.ndarray[DTYPE_t, ndim=2] p_xh_out = np.zeros( shape(p_xh_0), dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=1] mu_arr_out = np.zeros( shape(mu_arr_in), dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=1] sig_arr_out = np.zeros( shape(sig_arr_in), dtype=DTYPE )
    cdef DTYPE_t c_denomA
    cdef DTYPE_t c_numerA
    cdef DTYPE_t c_denom_p_xh
    cdef DTYPE_t c_denomB
    cdef DTYPE_t c_numer_mean
    cdef DTYPE_t c_numer_std
    cdef DTYPE_t pi_sum
    cdef int i, j, t
    #denominator_A = 0.0
  
  
    cdef np.ndarray[DTYPE_t, ndim=1] p_init_out = np.zeros( N, dtype=DTYPE )
    # pi, p_xh, mu, sig
    for i in range(N):
      #pi
      if gamma_arr[0,i] >= 0:
        p_init_out[i] = gamma_arr[0,i]
      else: # negative or gamma is nan
        print( ' pi[t=%d] is negative !' % i )
        print( ' gamma_arr[0,i] = %g' % gamma_arr[0,i] )
        p_init_out[i] = 0
  
      # p_xh
      c_denomA = 0.0
      #denominator_A = 0.0
      for t in range(T-1):
        c_denomA += gamma_arr[t, i]
  
      for j in range(N):
        c_numerA = 0.0
        for t in range(T-1):
          c_numerA += xi_arr[t, i, j]
  
        if c_denomA == 0.0 :
          """
          c_numerA is also zero. 
          So, it is safe to set c_denomA = 1.0
          For zero transition case, asstign p_xh_out[i,i] = 1.0
          """
          if i == j:
            p_xh_out[i,j] = 1.0 
          else:
            p_xh_out[i,j] = 0.0
        else:
          """
          Usual
          """
          p_xh_out[i, j] = c_numerA / c_denomA
  
        # Make it non ngative
        if p_xh_out[i, j] < 0:
          print( 'p_xh_out[%d, %d] is negative!!' % (i, j) )
          p_xh_out[i, j] = 0
  
      # p_xh row normalization 
      #denominator_p_xh = sum( p_xh_out[i,:] )
      c_denom_p_xh = sum( p_xh_out[i,:] )
      #if denominator_p_xh != 0:
      if c_denom_p_xh != 0:
        p_xh_out[i,:] /= c_denom_p_xh # unnecessary step... actually
      else:
        p_xh_out[i,:] = np.ones( shape(p_xh_out[i,:], dtype=DTYPE ) / (1.0 * N) )
  
      # mu
      c_denomB= c_denomA + gamma_arr[T-1, i]
      #denominator_B = denominator_A + gamma_arr[T-1, i]
  
      c_numer_mean = 0.0
      c_numer_std = 0.0
      for t in range(T):
        c_numer_mean += gamma_arr[t, i] * o_arr[t]
        c_numer_std  += gamma_arr[t, i] \
                           * ( o_arr[t] - mu_arr_in[i] ) \
                           * ( o_arr[t] - mu_arr_in[i] )
  
      mu_arr_out[i]  = c_numer_mean / c_denomB
      sig_arr_out[i] = sqrt( c_numer_std / c_denomB )
      if sig_arr_out[i] < sig_min:
        sig_arr_out[i] = sig_min
  
  
    # pi normalization
    pi_sum=0
    for i in range(N):
      pi_sum += p_init_out[i]

    if pi_sum != 0:
      p_init_out /= pi_sum
    else:
      p_init_out = 1.0/N * np.ones( shape( p_init_out ), dtype=DTYPE )
  
  
    return p_init_out, p_xh_out, mu_arr_out, sig_arr_out
  
  
  def baum_welch_1(np.ndarray[DTYPE_t, ndim=1] o_arr, 
                   np.ndarray[DTYPE_t, ndim=1] p_init_in, 
                   np.ndarray[DTYPE_t, ndim=2] p_xh_in, 
                   np.ndarray[DTYPE_t, ndim=1] mu_arr_in, 
                   np.ndarray[DTYPE_t, ndim=1] sig_arr_in): 
    alpha_arr_r, a_scale, log_probability = construct_scaled_alpha(
                                              o_arr,
                                              p_init_in,
                                              p_xh_in,
                                              mu_arr_in,
                                              sig_arr_in)
    beta_arr_r = construct_scaled_beta( 
                    o_arr,
                    p_init_in,
                    p_xh_in,
                    mu_arr_in,
                    sig_arr_in,
                    a_scale)
    gamma_arr = construct_gamma(alpha_arr_r, beta_arr_r)
    xi_arr = construct_xi( alpha_arr_r, 
                           beta_arr_r,
                           o_arr,
                           p_xh_in,
                           mu_arr_in,
                           sig_arr_in)
  
    p_init_out, p_xh_out, mu_arr_out, sig_arr_out = construct_params(
                                                      gamma_arr,
                                                      xi_arr,
                                                      o_arr,
                                                      mu_arr_in,
                                                      sig_arr_in)
  
    return p_init_out, p_xh_out, mu_arr_out, sig_arr_out, log_probability, \
           alpha_arr_r, beta_arr_r, gamma_arr, xi_arr
                              
    
  def hmm_param_estimation(
      np.ndarray[DTYPE_t, ndim=1] o_arr, 
      np.ndarray[DTYPE_t, ndim=1] p_init, 
      np.ndarray[DTYPE_t, ndim=2] p_xh, 
      np.ndarray[DTYPE_t, ndim=1] mu_arr, 
      np.ndarray[DTYPE_t, ndim=1] sig_arr,
      log_likelihood_cutoff=1e-4,
      max_iter=500,
      fix_emissions=False):
    """
    Execute baum_welch_1 recusively until meets the cutoff value.
    """
    cdef DTYPE_t log_prob1 = -99999.0
    cdef DTYPE_t log_prob2 = -9999999.0
    cdef int n_iter = 0
    while (fabs( log_prob1 - log_prob2 ) > log_likelihood_cutoff) \
        and (n_iter < max_iter):
      log_prob1 = log_prob2 # preveious log2 value assigned to log1
      p_init, p_xh, mu_arr, sig_arr, log_prob2, alpha_arr_r, beta_arr_r, gamma_arr, xi_arr = baum_welch_1(o_arr,
             p_init,
             p_xh,
             mu_arr,
             sig_arr)
      n_iter += 1
      print('%d th try, log_prob = %g' % (n_iter, log_prob2) )
      #print('%g' % (log_prob2 - log_prob1) )
  
    return p_init, p_xh, mu_arr, sig_arr, log_prob2, alpha_arr_r, beta_arr_r, gamma_arr, xi_arr, n_iter
   
  
  
  def viterbi(np.ndarray[DTYPE_t, ndim=1] o_arr, 
              np.ndarray[DTYPE_t, ndim=1] p_init, 
              np.ndarray[DTYPE_t, ndim=2] p_xh, 
              np.ndarray[DTYPE_t, ndim=1] mu_arr, 
              np.ndarray[DTYPE_t, ndim=1] sig_arr):
    """
    log-scaled viterbi
    ref) Sagemath, chmm.pyx
    """
    # init
    cdef np.ndarray[DTYPE_t, ndim=2] log_delta_arr = np.zeros( [T, N], dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=2] psi = np.zeros( [T, N], dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=2] log_p_xh = np.zeros( shape(p_xh), dtype=DTYPE )
    cdef np.ndarray[DTYPE_t, ndim=1] xh_arr_post = np.zeros( shape(o_arr), dtype=DTYPE )
    cdef int i, j, t
    cdef DTYPE_t mx, tmp, koko_max
    cdef int idx, koko_idx
    for i in range(N):
      for j in range(N):
        log_p_xh[i,j] = log( p_xh[i,j] )
    #log_p_xh = log(p_xh)
    for i in range(N):
      log_delta_arr[0, i] = log( p_init[i] ) \
                             + log( emmission_prob( o_arr[0], i, mu_arr[i], sig_arr[i] ) ) 
  
  
    # induction
    for t in range(1, T):
      for j in range(N):
        mx = -999999999999
        idx = -1
        for i in range(N):
          tmp = log_delta_arr[t-1, i] + log_p_xh[i, j]
          if tmp > mx:
            mx = tmp
            idx = i
        log_delta_arr[t, j] = mx + log( emmission_prob(o_arr[t], j, mu_arr[j], sig_arr[j]) )
        psi[t, j] = idx
  
    #termination
    koko_max = -999999
    koko_idx = -1
    for i in range(N):
      if log_delta_arr[T-1, i] > koko_max:
        koko_max = log_delta_arr[T-1, i]
        koko_idx = i
  
    # Back tracking
    xh_arr_post[T-1] = DTYPE (koko_idx)
    for t in range(T-2, -1, -1): # T-2, T-1, ..., 0
      xh_arr_post[t] = psi[t+1, xh_arr_post[t+1]]
  
    return xh_arr_post, mx, log_delta_arr, psi
  
  def construct_x_arr(np.ndarray[DTYPE_t, ndim=1] xh_arr_post, 
                      np.ndarray[DTYPE_t, ndim=1] mu_arr_post):
    cdef int t
    cdef np.ndarray[DTYPE_t, ndim=1] x_arr = np.zeros( shape(xh_arr_post), dtype=DTYPE )
    for t in range(T):
      x_arr[t] = mu_arr_post[xh_arr_post[t]]
  
    return x_arr
  
  # Main One time param estimation
  p_init_out, p_xh_out, mu_arr_out, sig_arr_out, log_probability, alpha_arr_r, beta_arr_r, gamma_arr, xi_arr = baum_welch_1(o_arr,
                       p_init_0,
                       p_xh_0,
                       mu_arr_0,
                       sig_arr_0)
  
  #p_init_out, p_xh_out, mu_arr_out, sig_arr_out, log_probability, alpha_arr_r, beta_arr_r, gamma_arr, xi_arr = baum_welch_1(o_arr,
  #           p_init_out,
  #           p_xh_out,
  #           mu_arr_out,
  #           sig_arr_out)
  #
  ### Parameter estimation
  #show()
  p_init_out, p_xh_out, mu_arr_out, sig_arr_out, log_probability, alpha_arr_r, beta_arr_r, gamma_arr, xi_arr, n_iter = hmm_param_estimation(o_arr,
             p_init_out,
             p_xh_out,
             mu_arr_out,
             sig_arr_out)
  
  
  xh_arr_post, mx, log_delta_arr, psi = viterbi(o_arr, p_init_out, p_xh_out, mu_arr_out, sig_arr_out)
  
  x_arr = construct_x_arr(xh_arr_post, mu_arr_out)

  return log_probability, n_iter+1, xh_arr_post, x_arr, p_init_out, p_xh_out, mu_arr_out, sig_arr_out # +1 due to initial baum-welch
