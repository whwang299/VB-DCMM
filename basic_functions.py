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
## 0.

from os import path
from copy import deepcopy
from matplotlib.pylab import *
from scipy.optimize import curve_fit
from scipy.optimize import fmin_l_bfgs_b
from c_hmm_func import c_hmm_func
import random

from scipy.special import psi # Digamma func
import scipy.special as sc_special 
gamma_func = sc_special.gamma
from scipy.special import gammaln
# To avoid overloading with numpy.random.gamma
 
 
def kl_dirichlet( u1_arr, u2_arr ):
  """
  From
  eq. 68
  Ji et al., Variational Bayes for Continous Hidden Markov...

  u1_arr: param for first dirchlet dist.
  u2_arr: param for second dirchlet dist.
  """
  # If lists come, convert them to numpy array
  u1_arr = array( u1_arr )
  u2_arr = array( u2_arr )

  u10 = sum( u1_arr )
  u20 = sum( u2_arr )

  #s1 = log( gamma_func(u10) / gamma_func(u20) ) # overflow..
  #s2 = sum( log( gamma_func(u2_arr) / gamma_func(u1_arr) ) )
  s1 = gammaln( u10 ) - gammaln( u20 ) 
  s2 = sum( gammaln( u2_arr ) - gammaln( u1_arr ) )
  s3 = sum( (u1_arr - u2_arr)*( psi(u1_arr) - psi(u10) ) )

  return s1 + s2 + s3


def normalize_transition_matrix( A ):
  """
  A: transition matrix whose row-sum is 1
  """
  nA = deepcopy( A )
  K1, K2 = shape( A )
  if K1 != K2:
    raise Exception('A should be squre matrx!')
  for i in range( K1 ):
    row_vec = A[i,:]
    tot = sum(row_vec)
    nA[i,:] = row_vec/tot

  return nA

def normalize_list_of_transition_matrix( A_list ):
  """
  A_list = [A1, A2, A3, ..., ]
  """
  L = len( A_list )
  nA_list = deepcopy( A_list )
  for l in range(L):
    nA_list[l] = normalize_transition_matrix( A_list[l] )

  return nA_list





