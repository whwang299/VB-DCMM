# 2015 Feb 27.
# Wonseok Hwang
# Code written by Wonseok Hwang 
# License: GPLv3
# Given one data, do the vbDCMM analysis
# All idx starts from 0. 

# 0. Init.
import os
from matplotlib.pylab import *
from numpy.random import dirichlet
from copy import deepcopy
from vbDCMM_random_1_func import vbDCMM_random_1_func

def vbDCMM_global_func(
                xh_arr_post, p_xh_post, mu_arr_post, # From pre-hmm analysis
                K_min=1, 
                K_max=3, 
                uphi_prior=1.0,
                ua_prior=1.0, uad_prior=10.0,
                ub_prior=1.0, ubd_prior=10.0,
                ua_prior_ic=1.0, uad_prior_ic=20.0,
                ub_prior_ic=1.0, ubd_prior_ic=20.0,
                e_cutoff=1e-3, 
                max_iter_vb=500, 
                answer_fix_init_param='n',
                N_gl=10, # Trial number for global maximum searching
                answer_ub_from_hmm='n'): 
  """
  Inputs:

  xh_arr_post: Sequence of noise-filtered observable state.
  p_xh_post: Estimated transition matrix for observable states using HMM.
  mu_arr_post: Estimated average value of each observable state.
  uphi_prior, ua_prior, uad_prior, ub_prior, ubd_prior: prior parameters for initial parameter distribution. See our paper for the details.
  u**_prior_ic: prior parameters used in generation of random initial conditions.
  N_gl: Number of repeat (to avoid local minimum)
  answer_ub_from_hmm: if 'y', p_xh_post is used to make ub_prior and ubd_prior. See our paper for the details.



  Outputs:

  xhh_arr_post_best_list: Most probable sequence of internal state for each model.
  phi_star_arr_best_list: Estimated phi.
  a_star_arr_best_list: Estimated A (transition matrix for internal state) for each model.
  b_star_arr_best_list: Estimated B (transitoin matrix for observable state) for each model.
  Fphi_arr_best_list: Fphi for each model.
  Fa_arr_best_list: Fa for each model.
  Fb_arr_best_list: Fb for each model.
  logL_arr_best_list: Likelihood for each model.
  F_arr_best_list: F for each model.
  flag_coverged_best

  Following returns are same with above but contain the results of all try.

  xhh_arr_post_list_all_try
  phi_star_arr_list_all_try
  a_star_arr_list_all_try 
  b_star_arr_list_all_try
  Fphi_arr_list_all_try
  Fa_arr_list_all_try
  Fb_arr_list_all_try
  logL_arr_list_all_try
  F_arr_list_all_try
  flag_coverged_list_all_try

  """
  # 1.1. Varying options
  T = len(xh_arr_post)
  L = shape( p_xh_post )[0]
  # to include smart initial guess
  if L == 2:
    # Aligning according to mu
    if mu_arr_post[0] < mu_arr_post[1]:
      k1 = p_xh_post[0][1] # k1: low fret --> high fret
      k1i = p_xh_post[1][0] # k1i: high_fret --> low fret
    else:
      k1 = p_xh_post[1][0]
      k1i = p_xh_post[0][1]

  # Construct time list "t_set"(L x L matrix) for OmegaB
  # template
  t_set = []
  for i in range(L):
    t_set.append( [] )
    for j in range(L):
      t_set[i].append( [] )
  
  # t_set construction
  for t in range(T-1):
    i1 = int( xh_arr_post[t]   )
    i2 = int( xh_arr_post[t+1] )
  
    t_set[i1][i2].append(t)
  
  idx = 0
  Fphi_arr_best_list = []
  Fa_arr_best_list = []
  Fb_arr_best_list = []
  logL_arr_best_list = []
  F_arr_best_list = []

  xhh_arr_post_best_list = []
  phi_star_arr_best_list = []
  a_star_arr_best_list = []
  b_star_arr_best_list = []

  Fphi_arr_list_all_try = []
  Fa_arr_list_all_try = []
  Fb_arr_list_all_try = []
  logL_arr_list_all_try = []
  F_arr_list_all_try = []

  flag_coverged_list_all_try = []

  xhh_arr_post_list_all_try = []
  phi_star_arr_list_all_try = []
  a_star_arr_list_all_try = []
  b_star_arr_list_all_try = []

  i_k = -1
  for K in range(K_min, K_max+1):
    i_k += 1

    Fphi_arr_list_all_try.append( [] )
    Fa_arr_list_all_try.append( [] )
    Fb_arr_list_all_try.append( [] )
    logL_arr_list_all_try.append( [] )
    F_arr_list_all_try.append( [] )

    flag_coverged_list_all_try.append( [] )
  
    xhh_arr_post_list_all_try.append( [] )
    phi_star_arr_list_all_try.append( [] )
    a_star_arr_list_all_try.append( [] )
    b_star_arr_list_all_try.append( [] )

    # To find global maximum, try multiple # with random initial conidition
    Fend_arr = []
    for i_gl in range( N_gl ):
      xhh_arr_post, phi_star_arr, a_star_arr, b_star_arr, Fphi_arr, Fa_arr, Fb_arr, logL_arr, F_arr, flag_coverged = vbDCMM_random_1_func(
                  xh_arr_post, p_xh_post, mu_arr_post, # From pre-hmm analysis
                  t_set, # t_set from xh_arr_post
                  K=K, 
                  uphi_prior=uphi_prior,
                  ua_prior=ua_prior, uad_prior=uad_prior,
                  ub_prior=ub_prior, ubd_prior=ubd_prior,
                  ua_prior_ic=ua_prior_ic, uad_prior_ic=uad_prior_ic,
                  ub_prior_ic=ub_prior_ic, ubd_prior_ic=ubd_prior_ic,
                  e_cutoff=e_cutoff, 
                  max_iter_vb=max_iter_vb,
                  answer_ub_from_hmm=answer_ub_from_hmm)
         
      xhh_arr_post_list_all_try[i_k].append( xhh_arr_post )
      phi_star_arr_list_all_try[i_k].append(phi_star_arr)
      a_star_arr_list_all_try[i_k].append(a_star_arr)
      b_star_arr_list_all_try[i_k].append(b_star_arr)
  
      Fphi_arr_list_all_try[i_k].append( Fphi_arr )
      Fa_arr_list_all_try[i_k].append( Fa_arr ) 
      Fb_arr_list_all_try[i_k].append( Fb_arr ) 
      logL_arr_list_all_try[i_k].append( logL_arr)
      F_arr_list_all_try[i_k].append( F_arr )

      flag_coverged_list_all_try[i_k].append( flag_coverged )

      Fend_arr.append( F_arr[-1] )

    # Pick Global maximum
    _idx = argmax( Fend_arr )
    xhh_arr_post_best = xhh_arr_post_list_all_try[i_k][_idx]
    phi_star_arr_best = phi_star_arr_list_all_try[i_k][_idx]
    a_star_arr_best = a_star_arr_list_all_try[i_k][_idx]
    b_star_arr_best = b_star_arr_list_all_try[i_k][_idx]
  
    Fphi_arr_best = Fphi_arr_list_all_try[i_k][_idx]
    Fa_arr_best = Fa_arr_list_all_try[i_k][_idx]
    Fb_arr_best = Fb_arr_list_all_try[i_k][_idx]
    logL_arr_best = logL_arr_list_all_try[i_k][_idx]
    F_arr_best = F_arr_list_all_try[i_k][_idx]

    flag_coverged_best = flag_coverged_list_all_try[i_k][_idx]

    # Add best result
    xhh_arr_post_best_list.append( xhh_arr_post_best )
    phi_star_arr_best_list.append(phi_star_arr_best)
    a_star_arr_best_list.append(a_star_arr_best)
    b_star_arr_best_list.append(b_star_arr_best)

    Fphi_arr_best_list.append( Fphi_arr_best )
    Fa_arr_best_list.append( Fa_arr_best ) 
    Fb_arr_best_list.append( Fb_arr_best ) 
    logL_arr_best_list.append( logL_arr_best )
    F_arr_best_list.append( F_arr_best )

  return xhh_arr_post_best_list, \
         phi_star_arr_best_list, \
         a_star_arr_best_list, \
         b_star_arr_best_list, \
         Fphi_arr_best_list, \
         Fa_arr_best_list, \
         Fb_arr_best_list, \
         logL_arr_best_list, \
         F_arr_best_list, \
         flag_coverged_best, \
         xhh_arr_post_list_all_try, \
         phi_star_arr_list_all_try, \
         a_star_arr_list_all_try, \
         b_star_arr_list_all_try, \
         Fphi_arr_list_all_try, \
         Fa_arr_list_all_try, \
         Fb_arr_list_all_try, \
         logL_arr_list_all_try, \
         F_arr_list_all_try, \
         flag_coverged_list_all_try


 
