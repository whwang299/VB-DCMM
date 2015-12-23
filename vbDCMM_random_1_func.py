# 2015 Feb 27.
# Wonseok Hwang
# License: GPLv3 or Manner.
# Given single data, perform the vbDCMM analysis
# All idx starts from 0. 
# Test on multiple number of states

# 0. Init.
import os
from matplotlib.pylab import *
from scipy.special import psi # Digamma func
from basic_functions import kl_dirichlet
from numpy.random import dirichlet
from c_hmm_func import c_hmm_func
from c_dcmm4vb_func import c_dcmm4vb_func
from c_dcmm_viterbi_func import c_dcmm_viterbi_func
from copy import deepcopy
#close('all')


def vbDCMM_random_1_func( xh_arr_post, p_xh_post, mu_arr_post, # From pre-hmm analysis
                t_set, # t_set from xh_arr_post
                K=2, 
                uphi_prior=1.0,
                ua_prior=1.0, uad_prior=10.0,
                ub_prior=1.0, ubd_prior=10.0,
                ua_prior_ic=1.0, uad_prior_ic=20.0,
                ub_prior_ic=1.0, ubd_prior_ic=20.0,
                e_cutoff=1e-3, 
                max_iter_vb=500,
                answer_ub_from_hmm='n'):
  """
  Inputs:

  xh_arr_post: Sequence of noise-filtered observable states
  p_xh_post: Estimated transition matrix for observable states using HMM.
  mu_arr_post: Estimated average value of each observable state.
  t_set: For calculation of W_b. See the paper.
  uphi_prior, ua_prior, uad_prior, ub_prior, ubd_prior: prior parameters for initial parameter distribution. See our paper for the details
  u**_prior_ic: prior parameters used in generation of random initial conditions.
  e_cutoff: Target precision for F.
  answer_ub_from_hmm: if 'y', use p_xh_post is used to make ub_prior and ubd_prior. See our paper for the details.



  Outputs
  xhh_arr_post: Estimated sequence of internal state in the model with K internal states.
  phi_star_arr: Estimated phi
  a_star_arr: Estimated transition matrix A for internal state
  b_star_arr: Estimated transition matrix B for observable state
  Fphi_arr: F(phi)
  Fa_arr: F(A)
  Fb_arr: F(B)
  logL_arr: log(likelihood)
  F_arr: F
  flag_coverged
  """
  # 1.1. Varying options
  T = len(xh_arr_post)
  L = shape( p_xh_post )[0]
  if L == 2:
    # Align according to mu
    if mu_arr_post[0] < mu_arr_post[1]:
      k1 = p_xh_post[0][1] # k1: low fret --> high fret
      k1i = p_xh_post[1][0] # k1i: high_fret --> low fret
    else:
      k1 = p_xh_post[1][0]
      k1i = p_xh_post[0][1]

  phi = zeros( K )
  A = zeros( [K, K] )
  B = zeros( [K, L, L] ) 
  
  # <Scheme> #
  # 0. Initial guess on phi, a, b -> q(z) #
  # Main loop start #############################################
  # 1. Estimate wpi, wA, wB,                                    # 
  # 2. Update W = u + w. u is Prior.                            #
  # 3. Update phi_star, a_star, b_star.                         #
  # 4. Calculate lower bound F
  # 5. Return to 1.                                             # 
  # Main loop end ###############################################
  # 3.1. Init guess on phi, A, B: q(z).
  uphi_arr = uphi_prior * ones( K )
  ua_arr = (uad_prior - ua_prior) * eye( K ) + ua_prior * ones( [K, K] )
  ua_arr_ic = (uad_prior_ic - ua_prior_ic) * eye( K ) + ua_prior_ic * ones( [K, K] )
  ub_arr = zeros([K,L,L])
  ub_arr_ic = zeros([K,L,L])
  if answer_ub_from_hmm == 'y':
    # Construct normed p_xh_post
    for k in range(K):
      for l1 in range(L):
        rate_min = min( p_xh_post[l1,:] )
        ub_arr[k,l1,:] = p_xh_post[l1,:] / rate_min

      ub_arr_ic[k,:,:] = (ubd_prior_ic - ub_prior_ic) * eye( L ) + ub_prior_ic * ones( [L, L] )
    #print('ub_arr', ub_arr)
  else:
    for k in range(K):
      ub_arr[k,:,:] = (ubd_prior - ub_prior) * eye( L ) + ub_prior * ones( [L, L] )
      ub_arr_ic[k,:,:] = (ubd_prior_ic - ub_prior_ic) * eye( L ) + ub_prior_ic * ones( [L, L] )
 
  A = array( [  dirichlet( ua_arr_ic[0,:] )  ] )
  for k in range(K):
    # A init
    if k > 0: # First one already made
      A = append(A, array( [ dirichlet( ua_arr_ic[k,:] ) ] ), 0)

    # B init
    _B = array( [  dirichlet( ub_arr_ic[k,0,:] )  ] )
    for l in range(1, L):
      _B = append(_B, array( [ dirichlet( ub_arr_ic[k,l,:] ) ] ), 0)
    B[k,:,:] = _B 
    # phi init
    phi[k] = 1.0/K


  # 4. Main loop: q(phi) -> q(z) -> q(phi) ...
  # templates
  wphi_arr = zeros( [K] )
  wa_arr = zeros( [T-2, K, K] )
  wb_arr = zeros( [T-1, K] )
  
  OmegaB_arr = zeros( [K, L, L] )
  
  Wphi_arr = zeros( [K] )
  Wa_arr = zeros( [K, K] )
  Wb_arr = zeros( [K, L, L] )
  
  phi_star_arr = zeros( [K] )
  a_star_arr = zeros( [K, K] )
  b_star_arr = zeros( [K, L, L] ) # Change K position for consistency with c_dcmm

  gamma_arr, xi_arr, logL = c_dcmm4vb_func(xh_arr_post, phi, A, B)

  #print(gamma_arr, xi_arr, logL)
  gamma_arr = deepcopy( gamma_arr[1:, :] )
  xi_arr = deepcopy( xi_arr[1:T-1,:,:] )
  _min = -1e+10
  Fphi_arr = [_min]
  Fa_arr = [_min]
  Fb_arr = [_min]
  logL_arr = [_min]
  F_arr = [_min]
  flag_coverged = 0
  for i_vb in range(max_iter_vb):
    # 3.2. Estimate wpi, wA, wB
    # gamma_arr[0] = blank
    # xi_arr[0,:,:], xi_arr[T-1,:,:] = blank
    # In c_dcmm, z: t=1, ... t=T-1, idx:0..T-2
    # In c_dcmm, x: t=0, ... t=T-1, idx:0..T-1
    # Here, z: t=0, ... t=T-2, idx:0..T-2
    # Here, x: t=0, ... t=T-1, idx:0..T-1
    # redefine gamma_arr, xi_arr to make consistent idx with
    # this script 
    # 1. Estimate wpi, wA, wB,                                    # 
    wphi_arr[:] = gamma_arr[0,:]
    wa_arr[:] = xi_arr[:]
    wb_arr[:] = gamma_arr[:]
  
    # 1.2. OmegaB
    OmegaB_arr[:] = 0
    for l1 in range(L):
      for l2 in range(L):
        t_sub_set = t_set[l1][l2]
        for k in range(K): 
          for t in t_sub_set:
            OmegaB_arr[k, l1, l2] += wb_arr[t, k] 
    
    
    # 2. Update W = u + w. u is Prior.                            #
    Wphi_arr[:] = wphi_arr[:] + uphi_arr[:]
    Wa_arr[:] = sum( wa_arr, axis=0) + ua_arr[:] # sum along the time
    Wb_arr[:] = OmegaB_arr[:] + ub_arr[:]

    #print('Wb_arr', OmegaB_arr)
    # 3. Update phi_star, a_star, b_star.                         #
    phi_star_arr[:] = exp( psi( Wphi_arr ) - psi( sum( Wphi_arr ) ) )
    for k in range(K):
      a_star_arr[k,:] = exp( psi(Wa_arr[k,:]) - psi( sum(Wa_arr[k,:]) ) )
      for l in range(L):
        b_star_arr[k, l, :] = exp( psi( Wb_arr[k, l, :]) - psi( sum(Wb_arr[k, l, :]) ) )

    # Renormalize phi-star_arr
    phi_star_arr[:] = phi_star_arr[:] / sum( phi_star_arr[:] )
  
    # c_dcmmvb_func again
    gamma_arr, xi_arr, logL = c_dcmm4vb_func(xh_arr_post, phi_star_arr, a_star_arr, b_star_arr)

    gamma_arr = deepcopy( gamma_arr[1:, :] )
    xi_arr = deepcopy( xi_arr[1:T-1,:,:] )

       
    #print('Wphi: ', Wphi_arr, ' Wa: ', Wa_arr, ' Wb: ', Wb_arr)
    # 4. Calculate lower bound F
    Fphi = -kl_dirichlet( Wphi_arr, uphi_arr )
    Fa = 0
    for k in range( K ):
      Fa += -kl_dirichlet( Wa_arr[k,:], ua_arr[k,:] )

    Fb = 0
    for k in range(K):
      for l in range(L):
        Fb += -kl_dirichlet( Wb_arr[k, l, :], ub_arr[k, l, :] )

    F = Fphi + Fa + Fb + logL

    Fphi_arr.append( Fphi )
    Fa_arr.append( Fa )
    Fb_arr.append( Fb )
    logL_arr.append( logL )
    F_arr.append( F )

    _resi_F = F_arr[i_vb+1] - F_arr[i_vb]
    _resi_Fphi = Fphi_arr[i_vb+1] - Fphi_arr[i_vb]
    _resi_Fa = Fa_arr[i_vb+1] - Fa_arr[i_vb]
    _resi_Fb = Fb_arr[i_vb+1] - Fb_arr[i_vb]
    _resi_logL = logL_arr[i_vb+1] - logL_arr[i_vb]
    print('vbDCMM K=%d, %dth iter, F = %g, dF = %g' % (K, i_vb, F_arr[i_vb+1],  _resi_F ) )

    #print(_resi_Fa, _resi_Fb, _resi_Fphi, _resi_logL)
    if abs( _resi_F) < e_cutoff:
      if abs( _resi_logL) < e_cutoff:
        if abs(_resi_Fphi) < e_cutoff:
          if abs( _resi_Fa) < e_cutoff:
            if abs( _resi_Fb) < e_cutoff:
              flag_coverged = 1
              break

  # 5. Return to 1.                                             # 
  # Main loop end 
  
  # Viterbi for DCMM
  xhh_arr_post, log_p_o_xstar =  c_dcmm_viterbi_func( xh_arr_post, phi_star_arr, a_star_arr, b_star_arr ) 
  xhh_arr_post = xhh_arr_post[1:]

  # Trim initial artificial min.
  Fphi_arr = Fphi_arr[1:]
  Fa_arr = Fa_arr[1:]
  Fb_arr = Fb_arr[1:]
  logL_arr = logL_arr[1:]
  F_arr = F_arr[1:]
  # Multimodal factor shall be added later.

  return xhh_arr_post, phi_star_arr, a_star_arr, b_star_arr, Fphi_arr, Fa_arr, Fb_arr, logL_arr, F_arr, flag_coverged
 
