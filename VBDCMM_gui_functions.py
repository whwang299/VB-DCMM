# Wonseok Hwang
# GPLv3
from matplotlib.pylab import *
from copy import deepcopy
import tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
import scipy.special as sc_special 
gamma_func = sc_special.gamma
from scipy.misc import comb
from basic_functions import normalize_transition_matrix
from basic_functions import normalize_list_of_transition_matrix


from c_hmm_func import c_hmm_func
from vbDCMM_global_func import vbDCMM_global_func
from vbDCMM_random_1_func import vbDCMM_random_1_func

import matplotlib
#matplotlib.rcParams.update({'figure.autolayout': True})

def VBDCMM_gui_analaysis_func(o_arr):
  """
  """

  # class
  class VBDCMM_analysis:
    def __init__(self, parent, o_arr):
      # ---------- params ----------
      self.parent = parent
      self.o_arr = deepcopy( o_arr )
      self.T_o = len( o_arr )
      self.T_xhh = self.T_o - 1
      self.K_esti_max = 4
      self.L_esti_max = 4

      
      self.time_o = arange( 0, self.T_o, 1 )
      self.time_xhh = arange( 0, self.T_xhh, 1 )


      # ---------- Grids coordinates ----------
      grid_label_L_esti = (0, 0)  # row, column, rowspan, columnspan
      grid_spinbox_L_esti = (0, 1)
      grid_label_mu_hmm = (0, 2)
      grid_entry_mu1_hmm = (0, 3)
      grid_entry_mu2_hmm = (1, 3)
      grid_entry_mu3_hmm = (2, 3)
      grid_entry_mu4_hmm = (3, 3)
      grid_label_sig_hmm = (0, 4)
      grid_entry_sig1_hmm = (0, 5)
      grid_entry_sig2_hmm = (1, 5)
      grid_entry_sig3_hmm = (2, 5)
      grid_entry_sig4_hmm = (3, 5)
      grid_label_B_hmm = (0, 6)
      grid_entry_B_hmm = (0, 7)
      grid_label_mu_esti = (0, 7+self.L_esti_max)
      grid_entry_mu1_esti = (0, 7+self.L_esti_max+1)
      grid_entry_mu2_esti = (0+1, 7+self.L_esti_max+1)
      grid_entry_mu3_esti = (0+2, 7+self.L_esti_max+1)
      grid_entry_mu4_esti = (0+3, 7+self.L_esti_max+1)
      grid_label_sig_esti = (0, 7+self.K_esti_max+2)
      grid_entry_sig1_esti = (0, 7+self.K_esti_max+3)
      grid_entry_sig2_esti = (0+1, 7+self.K_esti_max+3)
      grid_entry_sig3_esti = (0+2, 7+self.K_esti_max+3)
      grid_entry_sig4_esti = (0+3, 7+self.K_esti_max+3)
      grid_label_B0_esti = (0, 7+self.K_esti_max+4)
      grid_entry_B0_esti = (0, 7+self.K_esti_max+5)


      grid_button_hmm_analysis = (0+self.L_esti_max, 2, 1, 2)
      grid_canvas_hmm = (0+self.L_esti_max, 0, 1, self.L_esti_max*4 + 6)


      grid_label_K_min= (0+self.L_esti_max+1, 0)  # row, column, rowspan, columnspan
      grid_spinbox_K_min= (0+self.L_esti_max+1, 1)
      grid_label_K_max= (0+self.L_esti_max+1, 2)  # row, column, rowspan, columnspan
      grid_spinbox_K_max= (0+self.L_esti_max+1, 3)
      grid_label_rough_gamma = (0+self.L_esti_max+1, 4) # Number of internal states transition?
      grid_entry_rough_gamma = (0+self.L_esti_max+1, 5) # Number of internal states transition?
      grid_label_n_repeats = (0+self.L_esti_max+1, 6)
      grid_entry_n_repeats = (0+self.L_esti_max+1, 7)
      grid_label_F_methods = (0+self.L_esti_max+1, 8) 
      grid_spinbox_F_methods = (0+self.L_esti_max+1, 9) 
      grid_button_VBDCMM_analysis = (0+self.L_esti_max+2, 2, 1, 2)
      grid_canvas_VBDCMM = (0+self.L_esti_max+2, 0, 1, self.L_esti_max*4+6)
      grid_canvas_F = (0+self.L_esti_max+2, 12, 1, self.L_esti_max*4+6)

      grid_label_K_esti = (0+self.L_esti_max+3, 0)
      grid_entry_K_esti = (0+self.L_esti_max+3, 1)
      grid_label_A_esti = (0+self.L_esti_max+3, 2) 
      grid_entry_A_esti = (0+self.L_esti_max+3, 3) 
      grid_label_B1_esti = (0+self.L_esti_max+3, 3+self.K_esti_max)
      grid_entry_B1_esti = (0+self.L_esti_max+3, 3+self.K_esti_max+1)
      grid_label_B2_esti = (0+self.L_esti_max+3, 3+self.K_esti_max+1+1*self.L_esti_max)
      grid_entry_B2_esti = (0+self.L_esti_max+3, 3+self.K_esti_max+2+1*self.L_esti_max)
      grid_label_B3_esti = (0+2*self.L_esti_max+3, 3+self.K_esti_max+0*self.L_esti_max)
      grid_entry_B3_esti = (0+2*self.L_esti_max+3, 3+self.K_esti_max+1+0*self.L_esti_max)
      grid_label_B4_esti = (0+2*self.L_esti_max+3, 3+self.K_esti_max+1+1*self.L_esti_max)
      grid_entry_B4_esti = (0+2*self.L_esti_max+3, 3+self.K_esti_max+2+1*self.L_esti_max)
      grid_label_save_filename = (0+self.L_esti_max+3+1, 3+self.K_esti_max+2+2*self.L_esti_max, 1,2 )
      grid_entry_save_filename = (0+self.L_esti_max+3+2, 3+self.K_esti_max+2+2*self.L_esti_max, 1,2 )
      grid_button_save = (0+self.L_esti_max+3+3, 3+self.K_esti_max+2+2*self.L_esti_max, 2, 2)
 
      # ---------- Variables & Initial value ----------
      #now define fvariables
      self.L_default = 2
      self.K_default = 2
      self.L_esti = Tk.IntVar(master=parent)
      self.L_esti.trace('w', self.on_L_change )
      self.L_esti.set( self.L_default )
      self.K_min = 1
      self.K_max = 3

      str_exec = ''
      for i in range(1,self.L_esti_max+1): # initial parameters
        str_exec += '\n'+'self.mu'+str(i)+'_hmm = Tk.DoubleVar(master=parent)'
        str_exec += '\n'+'self.mu'+str(i)+'_hmm.trace("w", self.on_mu_hmm_change)'
        str_exec += '\n'+'self.sig'+str(i)+'_hmm = Tk.DoubleVar(master=parent)'
        str_exec += '\n'+'self.sig'+str(i)+'_hmm.trace("w", self.on_sig_hmm_change)'
      exec(str_exec)
  
      str_exec = ''
      for i in range(1,self.L_default+1):
        str_exec += '\n'+'self.mu'+str(i)+'_hmm.set(' + str( 0.1 + 0.8/(self.L_default-1)* (i-1) ) + ')'
      exec(str_exec)
   
      str_exec = ''
      for i in range(1,self.L_default+1):
        str_exec += '\n'+'self.sig'+str(i)+'_hmm.set('+ str( 0.15/self.L_default ) + ')'
      exec(str_exec)
  
      str_exec = ''
      for i in range(1,self.L_esti_max+1):
        for j in range(1,self.L_esti_max+1):
          str_exec += '\n'+'self.B' + str(i) + str(j) + '_hmm = Tk.DoubleVar(master=parent)'
          str_exec += '\n'+'self.B0' + str(i) + str(j) + '_esti = Tk.DoubleVar(master=parent)'
          #str_exec += '\n'+'self.B' + str(i) + str(j) + '_hmm.trace("w", self.on_B_hmm_change)'
      for i in range(1,self.L_default+1):
        for j in range(1,self.L_default+1):
          if i != j:
            str_exec += '\n'+'self.B' + str(i) + str(j) + '_hmm.set(' + str( 0.01 ) + ')'
          if i == j:
            str_exec += '\n'+'self.B' + str(i) + str(j) + '_hmm.set(' + str( 1 -  (self.L_default-1)*0.01) + ')'
      exec(str_exec)

      self.K_min = Tk.IntVar( master = parent )
      self.K_min.set(1)
      self.K_max = Tk.IntVar( master = parent )
      self.K_max.set(3)

      self.rough_gamma = Tk.DoubleVar(master=parent)
      self.rough_gamma.set( 0.001 )

      self.n_repeats = Tk.IntVar(master=parent)
      self.n_repeats.set( 1 )

      self.F_methods = Tk.StringVar(master=parent)
      self.F_arr = zeros( 4 )

      # Result of analysis part
      self.K_esti = Tk.IntVar(master=parent)
      self.K_esti.set(0)
      str_exec = ''
 
      str_exec = ''
      for i in range(1, self.K_esti_max+1):
        for j in range(1, self.K_esti_max+1):
          str_exec += '\n'+'self.A' + str(i) + str(j) + '_esti = Tk.DoubleVar(master=parent)'
          str_exec += '\n'+'self.A' + str(i) + str(j) + '_esti.set(0)'
      exec(str_exec)
 
      str_exec = ''
      for i in range(1, self.L_esti_max+1):
        str_exec += '\n'+'self.mu'+str(i)+'_esti = Tk.DoubleVar(master=parent)'
        str_exec += '\n'+'self.sig'+str(i)+'_esti = Tk.DoubleVar(master=parent)'

        str_exec += '\n'+'self.mu'+str(i)+'_esti.set(0)'
        str_exec += '\n'+'self.sig'+str(i)+'_esti.set(0)'
      exec(str_exec)
      
      for k in range(1, self.K_esti_max+1):
        for i in range(1,self.L_esti_max+1):
          for j in range(1,self.L_esti_max+1):
            str_exec += '\n'+'self.B' + str(k) + str(i) + str(j) + '_esti = Tk.DoubleVar(master=parent)'
            str_exec += '\n'+'self.B' + str(k) + str(i) + str(j) + '_esti.set(0)'
      exec(str_exec)

      self.save_filename = Tk.StringVar(master=parent, value='batch1')
        
      # ---------- GUI Elements ----------
      self.label_L_esti = Tk.Label( master=parent, width=2, text='N:' )
      self.label_L_esti.grid( row=grid_label_L_esti[0], column=grid_label_L_esti[1] )
      self.spinbox_L_esti = Tk.Spinbox( parent, width=2, value=(1,2,3,4), textvariable=self.L_esti )
      self.L_esti.set( self.L_default )
      self.spinbox_L_esti.grid( row=grid_spinbox_L_esti[0], column=grid_spinbox_L_esti[1] )
      # mu and sig_hmm
      self.label_mu_hmm = Tk.Label( master=parent, width=2, text='mu:')
      self.label_mu_hmm.grid( row=grid_label_mu_hmm[0], column=grid_label_mu_hmm[1] )
      self.label_sig_hmm = Tk.Label( master=parent, width=2, text='sig:')
      self.label_sig_hmm.grid( row=grid_label_sig_hmm[0], column=grid_label_sig_hmm[1] )
      _s = ''
      for i in range(1,self.L_esti_max+1):
        _s += '\n'+'self.entry_mu' + str(i) + '_hmm = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.mu'+str(i)+'_hmm, width=5)'
        _s += '\n'+'self.entry_sig' + str(i) + '_hmm = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.sig'+str(i)+'_hmm, width=5)'
        if i <= self.L_default:
          _s += '\n'+'self.entry_mu' + str(i) + '_hmm.config(relief=Tk.SOLID)' 
          _s += '\n'+'self.entry_sig' + str(i) + '_hmm.config(relief=Tk.SOLID)' 
      exec(_s)

      for i in range(1,self.L_esti_max+1):
        _row = eval( 'grid_entry_mu' + str(i) + '_hmm[0]')
        _column = eval( 'grid_entry_mu' + str(i) + '_hmm[1]')
        _s += '\n'+'self.entry_mu' + str(i) + '_hmm.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
      exec(_s)
  
      _s = ''
      for i in range(1,self.L_esti_max+1):
        _row = eval( 'grid_entry_sig' + str(i) + '_hmm[0]')
        _column = eval( 'grid_entry_sig' + str(i) + '_hmm[1]')
        _s += '\n'+'self.entry_sig' + str(i) + '_hmm.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
      exec(_s)

      self.label_B_hmm = Tk.Label( master=parent, width=2, text='B:')
      self.label_B_hmm.grid( row=grid_label_B_hmm[0], column=grid_label_B_hmm[1] )

      self.label_B0_esti = Tk.Label( master=parent, width=4, text='B_esti:')
      self.label_B0_esti.grid( row=grid_label_B0_esti[0], column=grid_label_B0_esti[1] )


      _s = ''
      for i in range(1,self.L_esti_max+1):
        for j in range(1,self.L_esti_max+1):
          # "exec" doesnt work inside of loop.
          _row = (i-1) + grid_entry_B_hmm[0]
          _column = (j-1) + grid_entry_B_hmm[1]
          _s += '\n'+'self.entry_B' + str(i) + str(j) + '_hmm = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.B'+str(i)+str(j)+'_hmm, width=5)'
          if i <= self.L_default and j <= self.L_default:
            _s += '\n'+'self.entry_B' + str(i) + str(j) + '_hmm.config(relief=Tk.SOLID)'
 
          _s += '\n'+'self.entry_B' + str(i) + str(j) + '_hmm.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    

          _row = (i-1) + grid_entry_B0_esti[0]
          _column = (j-1) + grid_entry_B0_esti[1]
          _s += '\n'+'self.entry_B0' + str(i) + str(j) + '_esti = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.B0'+str(i)+str(j)+'_esti, width=5)'
          if i <= self.L_default and j <= self.L_default:
            _s += '\n'+'self.entry_B0' + str(i) + str(j) + '_esti.config(relief=Tk.SOLID)'
 
          _s += '\n'+'self.entry_B0' + str(i) + str(j) + '_esti.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    



      exec(_s)
      # ---------- Connect matplotlib canvas with Tk  ----------
    
      #self.button_hmm_analysis = Tk.Button( master=parent, 
      self.button_hmm_analysis = Tk.Button(master=parent, text='HMM!', width=8, command=self.on_click_button_hmm_analysis, height=4)
      self.button_hmm_analysis.grid( row=grid_button_hmm_analysis[0], column=grid_button_hmm_analysis[1],
                                          rowspan=grid_button_hmm_analysis[2], columnspan=grid_button_hmm_analysis[3],
                                          sticky='EW')

  

      self.label_K_min = Tk.Label( master=parent, width=5, text='K min:')
      self.label_K_min.grid( row=grid_label_K_min[0], column=grid_label_K_min[1] )

      self.spinbox_K_min = Tk.Spinbox(parent, width=2,  value=(1, 2, 3, 4), textvariable=self.K_min )
      self.spinbox_K_min.grid( row=grid_spinbox_K_min[0], column=grid_spinbox_K_min[1])
      self.K_min.set( 1 )

      self.label_K_max = Tk.Label( master=parent, width=5, text='K max:')
      self.label_K_max.grid( row=grid_label_K_max[0], column=grid_label_K_max[1] )


      self.spinbox_K_max = Tk.Spinbox(parent, width=2, value=(1, 2, 3, 4), textvariable=self.K_max )
      self.spinbox_K_max.grid( row=grid_spinbox_K_max[0], column=grid_spinbox_K_max[1])
      self.K_max.set( self.K_default )

      self.label_rough_gamma = Tk.Label( master=parent, width=6, text='gamma:')
      self.label_rough_gamma.grid( row=grid_label_rough_gamma[0], column=grid_label_rough_gamma[1] )

      self.entry_rough_gamma = Tk.Entry( parent, textvariable=self.rough_gamma, width=6 )
      self.entry_rough_gamma.grid( row=grid_entry_rough_gamma[0], column=grid_entry_rough_gamma[1] )


      self.label_n_repeats = Tk.Label( master=parent, width=8, text='n_repeats:')
      self.label_n_repeats.grid( row=grid_label_n_repeats[0], column=grid_label_n_repeats[1] )

      self.entry_n_repeats = Tk.Entry( parent, textvariable=self.n_repeats, width=2)
      self.entry_n_repeats.grid( row=grid_entry_n_repeats[0], column=grid_entry_n_repeats[1])

      self.label_F_methods = Tk.Label( master=parent, width=8, text='F methods:')
      self.label_F_methods.grid( row=grid_label_F_methods[0], column=grid_label_F_methods[1] )

      self.spinbox_F_methods = Tk.Spinbox( parent, width=8, value=('Based on G', 'Based on F') )
      self.spinbox_F_methods.grid( row=grid_spinbox_F_methods[0], column=grid_spinbox_F_methods[1] )

      self.button_VBDCMM_analysis = Tk.Button(master=parent, text='VB-DCMM!', width=8, command=self.on_click_button_VBDCMM_analysis, height=4)
      self.button_VBDCMM_analysis.grid( row=grid_button_VBDCMM_analysis[0], column=grid_button_VBDCMM_analysis[1],
                                          rowspan=grid_button_VBDCMM_analysis[2], columnspan=grid_button_VBDCMM_analysis[3],
                                          sticky='EW')


      self.label_save_filename = Tk.Label( master=parent, width=9, text='Filename:')
      self.label_save_filename.grid( row=grid_label_save_filename[0], column=grid_label_save_filename[1], 
          rowspan=grid_entry_save_filename[2],columnspan=grid_entry_save_filename[3])

      self.entry_save_filename = Tk.Entry(master=parent, width=9, textvariable=self.save_filename )
      self.entry_save_filename.grid( row=grid_entry_save_filename[0], column=grid_entry_save_filename[1],
          rowspan=grid_entry_save_filename[2],columnspan=grid_entry_save_filename[3])


      self.button_save = Tk.Button(master=parent, text='Save', width=4, command=self.on_click_button_save, height=4)
      self.button_save.grid( row=grid_button_save[0], column=grid_button_save[1],
                                          rowspan=grid_button_save[2], columnspan=grid_button_save[3],
                                          sticky='EW')

      self.label_K_esti = Tk.Label( master=parent, width=5, text='Best K:')
      self.label_K_esti.grid( row=grid_label_K_esti[0], column=grid_label_K_esti[1] )


      self.entry_K_esti = Tk.Entry(parent, width=5, textvariable=self.K_esti )
      self.entry_K_esti.grid( row=grid_entry_K_esti[0], column=grid_entry_K_esti[1] )

      self.label_A_esti = Tk.Label( master=parent, width=2, text='A:')
      self.label_A_esti.grid( row=grid_label_A_esti[0], column=grid_label_A_esti[1] )

      _s = ''
      for i in range(1,self.K_esti_max+1):
        for j in range(1,self.K_esti_max+1):
          _s += '\n'+'self.entry_A' + str(i) + str(j) + '_esti = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.A'+str(i)+str(j)+'_esti, width=5)'
          _row = (i-1) + grid_entry_A_esti[0]
          _column = (j-1) + grid_entry_A_esti[1]
          _s += '\n'+'self.entry_A' + str(i) + str(j) + '_esti.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
  
          if i <= self.K_default and j <= self.K_default:
            _s += '\n'+'self.entry_A' + str(i) + str(j) + '_esti.config(relief=Tk.SOLID)'
      exec(_s)



      self.label_mu_esti = Tk.Label( master=parent, width=6, text='mu_esti:')
      self.label_mu_esti.grid( row=grid_label_mu_esti[0], column=grid_label_mu_esti[1] )


      for i in range(1,self.L_esti_max+1):
        _row = eval( 'grid_entry_mu' + str(i) + '_esti[0]')
        _column = eval( 'grid_entry_mu' + str(i) + '_esti[1]')
        _s += '\n'+'self.entry_mu' + str(i) + '_esti = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.mu'+str(i)+'_esti, width=5)'
        _s += '\n'+'self.entry_mu' + str(i) + '_esti.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
        if i <= self.L_default:
          _s += '\n'+'self.entry_mu'+str(i)+'_esti.config(relief=Tk.SOLID)'

      exec(_s)
 
      self.label_sig_esti = Tk.Label( master=parent, width=6, text='sig_esti:')
      self.label_sig_esti.grid( row=grid_label_sig_esti[0], column=grid_label_sig_esti[1] )
 
      _s = ''
      for i in range(1,self.L_esti_max+1):
        _row = eval( 'grid_entry_sig' + str(i) + '_esti[0]')
        _column = eval( 'grid_entry_sig' + str(i) + '_esti[1]')
        _s += '\n'+'self.entry_sig' + str(i) + '_esti = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.sig'+str(i)+'_esti, width=5)'
        _s += '\n'+'self.entry_sig' + str(i) + '_esti.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
        if i <= self.L_default:
          _s += '\n'+'self.entry_sig'+str(i)+'_esti.config(relief=Tk.SOLID)'
      exec(_s)
  

      self.label_B1_esti = Tk.Label( master=parent, width=2, text='B1:')
      self.label_B2_esti = Tk.Label( master=parent, width=2, text='B2:')
      self.label_B3_esti = Tk.Label( master=parent, width=2, text='B3:')
      self.label_B4_esti = Tk.Label( master=parent, width=2, text='B4:')

      self.label_B1_esti.grid( row=grid_label_B1_esti[0], column=grid_label_B1_esti[1] )
      self.label_B2_esti.grid( row=grid_label_B2_esti[0], column=grid_label_B2_esti[1] )
      self.label_B3_esti.grid( row=grid_label_B3_esti[0], column=grid_label_B3_esti[1] )
      self.label_B4_esti.grid( row=grid_label_B4_esti[0], column=grid_label_B4_esti[1] )

      _s = ''
      for k in range(1, self.K_esti_max+1):
        for i in range(1,self.L_esti_max+1):
          for j in range(1,self.L_esti_max+1):
            # "exec" doesnt work inside of loop.
            _row = (i-1) + eval('grid_entry_B'+str(k)+'_esti[0]' )
            _column = (j-1) + eval('grid_entry_B'+str(k)+'_esti[1]' )
            _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '_esti = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.B'+str(k)+str(i)+str(j)+'_esti, width=5)'
            if k <= self.K_default and i <= self.L_default and j <= self.L_default:
              _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '_esti.config(relief=Tk.SOLID)'
 
            _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '_esti.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
      exec(_s)
      # ---------- Connect matplotlib canvas with Tk  ----------
      # -------- Canvas --------
      _figsize_x = 7
      _figsize_x_2 = 3
      _figsize_y = 2
      _padx_1 = 0
      _pady_1 = 0
      _rect = [0.10, 0.3, 0.8, 0.6]
      _rect_2 = [0.35, 0.3, 0.6, 0.55]
      fign_hmm = 101
      fign_VBDCMM = 102
      fign_F = 103

      self.x_arr_hmm = zeros( self.T_o )
      self.xhh_arr_VBDCMM = zeros( self.T_xhh )
      self.F_idx_arr = arange( self.K_min.get(), self.K_max.get()+1, 1)
      self.F_arr = zeros( self.K_max.get() - self.K_min.get() + 1 ) # used for either F or G
      self.Fend_arr = zeros( self.K_max.get() - self.K_min.get() + 1 )
      self.Gend_arr = zeros( self.K_max.get() - self.K_min.get() + 1 )

  
      self.fig_hmm = figure(num=fign_hmm, figsize=(_figsize_x, _figsize_y), facecolor='w')
      self.ax_hmm = self.fig_hmm.add_axes(_rect)
      self.ax_hmm.hold(False)
      self.ax_hmm.plot(self.time_o, self.o_arr, color=[0.7, 0.7, 0.7], linewidth=1)
      self.ax_hmm.set_ylim([-0.2, 1.2] )
      self.ax_hmm.set_yticks( [0, 0.5, 1] )
      self.ax_hmm.grid(True)
      self.ax_hmm.hold(False)
  
      self.canvas_hmm = FigureCanvasTkAgg(self.fig_hmm, master=parent)
      self.canvas_hmm.show()
      self.tk_canvas_hmm = self.canvas_hmm.get_tk_widget()
      self.tk_canvas_hmm.grid( row=grid_canvas_hmm[0], column=grid_canvas_hmm[1], 
                               rowspan=grid_canvas_hmm[2], columnspan=grid_canvas_hmm[3], 
                               padx=_padx_1, pady=_pady_1 )
      self.ax_hmm.set_xlabel('Time (dt=1)')
      self.ax_hmm.set_ylabel('O')
 
      cid1 = self.canvas_hmm.mpl_connect('button_press_event', self.on_figure_click)
      cid2 = self.canvas_hmm.mpl_connect('close_event', self.on_destroy)
   
      self.fig_VBDCMM = figure( num=fign_VBDCMM, figsize=(_figsize_x, _figsize_y), facecolor='w' )
      self.ax_VBDCMM = self.fig_VBDCMM.add_axes(_rect)
      self.ax_VBDCMM.hold(False)
      self.ax_VBDCMM.plot(self.time_xhh, self.xhh_arr_VBDCMM+1,'k', linewidth=2)
      self.ax_VBDCMM.set_ylim([0.5, self.K_default + 0.5] )
      self.ax_VBDCMM.set_yticks( arange(1, self.K_default+1, 1) )
      self.ax_VBDCMM.grid(True)
      self.ax_VBDCMM.hold(False)
 
      self.canvas_VBDCMM = FigureCanvasTkAgg(self.fig_VBDCMM, master=parent)
      self.canvas_VBDCMM.show()
      self.tk_canvas_VBDCMM = self.canvas_VBDCMM.get_tk_widget()
      self.tk_canvas_VBDCMM.grid( row=grid_canvas_VBDCMM[0], column=grid_canvas_VBDCMM[1], 
                             rowspan=grid_canvas_VBDCMM[2], columnspan=grid_canvas_VBDCMM[3], 
                             padx=_padx_1, pady=_pady_1 )
      self.ax_VBDCMM.set_xlabel('Time (dt=1)')
      self.ax_VBDCMM.set_ylabel('X')
 
      cid1 = self.canvas_VBDCMM.mpl_connect('button_press_event', self.on_figure_click)
      cid2 = self.canvas_VBDCMM.mpl_connect('close_event', self.on_destroy)
  
      self.fig_F = figure( num=fign_F, figsize=(_figsize_x_2, _figsize_y), facecolor='w' )
      self.ax_F = self.fig_F.add_axes(_rect_2)
      self.ax_F.plot(self.F_idx_arr, self.F_arr, '--rs')
      self.ax_F.grid(True)
      self.ax_F.hold(False)

      K_min = self.K_min.get()
      K_max = self.K_max.get()
      _xticks = arange( K_min, K_max+1, 1 )
      self.ax_F.set_xticks( _xticks )
      self.ax_F.set_xlim([K_min-0.5, K_max+0.5] )

      _yy = self.ax_F.get_ylim()
      _dyy = diff( _yy )
      _yy_min = _yy[0] - _dyy/7
      _yy_max = _yy[1] + _dyy/7
      self.ax_F.set_ylim([_yy_min, _yy_max] )
      self.ax_F.set_yticks(_yy)
      self.ax_F.set_xlabel('K')
      self.ax_F.set_ylabel('Evidence')


      self.canvas_F = FigureCanvasTkAgg(self.fig_F, master=parent)
      self.canvas_F.show()
      self.tk_canvas_F = self.canvas_F.get_tk_widget()
      self.tk_canvas_F.grid( row=grid_canvas_F[0], column=grid_canvas_F[1], 
                             rowspan=grid_canvas_F[2], columnspan=grid_canvas_F[3], 
                             padx=_padx_1, pady=_pady_1 )
  
      cid1 = self.canvas_F.mpl_connect('button_press_event', self.on_figure_click)
      cid2 = self.canvas_F.mpl_connect('close_event', self.on_destroy)
  
      # ---------- Callback functions ----------

    def on_destroy(self, event):
      """ 
      Without this..., loop after close the figure
      """
      print('Destroy!!') 
      self.parent.quit()
  
    def on_figure_click(self, event):
      print('you pressed', event.button, event.x, event.y, event.xdata, event.ydata)
      #print('e1_number: ' + str(self.e1_number.get()) )
  
    def on_click_button_save(self):
      print('Save traces %s.dat' % self.save_filename.get() )
      filename_0 = self.save_filename.get()

      # Save O_filt, 
      savetxt(filename_0 + '_O_filt_hmm.dat', self.x_arr_hmm, fmt='%.3f' )
      savetxt(filename_0 + '_X_best_VBDCMM.dat', self.xhh_arr_VBDCMM+1, fmt='%d' )
      
      savetxt(filename_0 + '_F_VBDCMM.dat', self.Fend_arr, fmt='%.3f' )
      savetxt(filename_0 + '_G_VBDCMM.dat', self.Gend_arr, fmt='%.3f' )
      
      # save A
      _K = self.K_esti.get()
      _T_xhh = self.T_xhh
      _T_o = self.T_o
      _L = self.L_esti.get()

      self._A = zeros( [_K, _K] ) 
      for i in range(1, _K+1):
        for j in range(1, _K+1):
          self._A[i-1, j-1] = eval( 'self.A' + str(i) + str(j) + '_esti.get()' )
      savetxt(filename_0+'_A_VBDCMM.dat', self._A, fmt='%.9f' )
      print('A')
      print(self._A)
  
      self._B = zeros( [_K, _L, _L] ) 
  
      # save B
      for k in range(1, _K+1):
        for i in range(1, _L+1):
          for j in range(1, _L+1):
            self._B[k-1, i-1, j-1] = eval( 'self.B' + str(k) + str(i) + str(j) + '_esti.get()' )
  
        savetxt(filename_0+'_B'+str(k)+'_VBDCMM.dat', self._B[k-1,:,:], fmt='%.9f' )
  
      # save mu and sig
      self._mu_arr = zeros( _L )
      self._sig_arr = zeros( _L )
  
      for i in range(1, _L+1):
        self._mu_arr[i-1] = eval('self.mu'+str(i)+'_esti.get()')
        self._sig_arr[i-1] = eval('self.sig'+str(i)+'_esti.get()')
  
      savetxt(filename_0+'_mu_arr_HMM.dat', self._mu_arr, fmt='%.3f' )
      savetxt(filename_0+'_sig_arr_HMM.dat', self._sig_arr, fmt='%.3f' )
  
    def on_click_button_hmm_analysis(self):
      print('HMM analysis')
      # Construct A matrix
      _L = self.L_esti.get()
      _T_xhh = self.T_xhh
      _T_o = self.T_o
      p_init_0_hmm = zeros( [_L, _L] )

      p_init_0_hmm = ones( _L ) / _L
      p_xh_0_hmm = zeros( [_L, _L] ) 
      mu_arr_0_hmm = zeros( _L )
      sig_arr_0_hmm = zeros( _L )
   
      for i in range(1, _L+1):
        mu_arr_0_hmm[i-1] = eval('self.mu'+str(i)+'_hmm.get()')
        sig_arr_0_hmm[i-1] = eval('self.sig'+str(i)+'_hmm.get()')

        for j in range(1, _L+1):
          p_xh_0_hmm[i-1, j-1] = eval( 'self.B' + str(i) + str(j) + '_hmm.get()' )
          
      log_probability_hmm, n_iter_hmm, xh_arr_post, x_arr, p_init_post_hmm, p_xh_post, mu_arr_post, sig_arr_post = c_hmm_func(
          self.o_arr,
          p_init_0_hmm,
          p_xh_0_hmm,
          mu_arr_0_hmm,
          sig_arr_0_hmm)

      self.x_arr_hmm = x_arr
      self.xh_arr_hmm = xh_arr_post
      self.B_hmm_a = normalize_transition_matrix( p_xh_post )# _a : after
      self.mu_arr_hmm_a = mu_arr_post
      self.sig_arr_hmm_a = sig_arr_post
      # Add some function for showing hmm analysis results
  
      self.ax_hmm.hold(False)
      self.ax_hmm.plot(self.time_o, self.o_arr, color=[0.7, 0.7, 0.7], linewidth=1)
      self.ax_hmm.hold(True)
      self.ax_hmm.plot(self.time_o, self.x_arr_hmm, color=[0, 0, 1], linewidth=1)
      
      self.ax_hmm.set_ylim([-0.2, 1.2] )
      self.ax_hmm.set_yticks( [0, 0.5, 1] )
      self.ax_hmm.grid(True)
      self.ax_hmm.set_xlabel('Time (dt=1)')
      self.ax_hmm.set_ylabel('O')
 
 
      self.fig_hmm.canvas.draw()

      # Update variables
      _s = ''
      for i in range(1,self.L_esti_max+1):
        if i <= _L:
          _s += '\n'+'self.mu' + str(i) + '_esti.set(' + str( self.mu_arr_hmm_a[i-1] ) +')'
          _s += '\n'+'self.sig' + str(i) + '_esti.set(' + str( self.sig_arr_hmm_a[i-1] ) +')'
          _s += '\n'+'self.entry_mu' + str(i) + '_esti.config(relief=Tk.SOLID)' 
          _s += '\n'+'self.entry_sig' + str(i) + '_esti.config(relief=Tk.SOLID)' 
        else:
          _s += '\n'+'self.mu' + str(i) + '_esti.set(0)' 
          _s += '\n'+'self.sig' + str(i) + '_esti.set(0)'
          _s += '\n'+'self.entry_mu' + str(i) + '_esti.config(relief=Tk.RIDGE)' 
          _s += '\n'+'self.entry_sig' + str(i) + '_esti.config(relief=Tk.RIDGE)' 
        for j in range(1, self.L_esti_max+1):
          if i <= _L and j <= _L:
            _s += '\n'+'self.B0' + str(i) + str(j) + '_esti.set(' + str( self.B_hmm_a[i-1,j-1] ) + ')'
            _s += '\n'+'self.entry_B0'  + str(i) + str(j) + '_esti.config(relief=Tk.SOLID)'
          else:
            _s += '\n'+'self.B0' + str(i) + str(j) + '_esti.set(' + str( 0 ) + ')'
            _s += '\n'+'self.entry_B0' + str(i) + str(j) + '_esti.config(relief=Tk.RIDGE)'

      exec(_s)
 
 
    def on_click_button_VBDCMM_analysis(self):
      print('VBDCMM analysis!!')
      # Construct A matrix
      _L = self.L_esti.get()
      _T_xhh = self.T_xhh
      _T_o = self.T_o

      xh_arr_post = self.xh_arr_hmm
      p_xh_post = self.B_hmm_a
      print(p_xh_post)
      mu_arr_post = self.mu_arr_hmm_a
      sig_arr_post = self.sig_arr_hmm_a
      K_min = self.K_min.get()
      K_max = self.K_max.get()
      self.F_idx_arr = arange( K_min, K_max+1, 1)
      uphi_prior = 1 
      ua_prior = 1
      uad_prior = 1/self.rough_gamma.get()
      print(uad_prior, uad_prior)
      ua_prior_ic = 0.3 # generate init. A dist
      uad_prior_ic = 200 # generate init. A dist
      ub_prior = 1 # not used 
      ubd_prior = 20 # not used 
      ub_prior_ic = 1 # used only to generate initial B distriubtion
      ubd_prior_ic = 20 # used only to generate initial B distribution

      answer_ub_from_hmm='y' # use hmm results for prior
      max_iter_vb = 300
      N_gl = self.n_repeats.get()
      k_factor_list = [  [[0.5,0.5],[1.5,1.5],[0.5,1.5],[1.5,0.5]], \
                     [[0.5,1.5],[1.5,0.5],[0.5,0.5],[1.5,1.5]], \
                     [[0.2,0.2],[1.2,1.2],[0.2,1.2],[1.2,0.2]], \
                     [[0.2,1.2],[1.2,0.2],[0.2,0.2],[1.2,1.2]], \
                     [[0.2,0.2],[2.0,2.0],[0.2,2.0],[2.0,0.2]], \
                     [[0.3,2.0],[2.0,0.3],[0.3,0.3],[2.0,2.0]]  ] #Kmax<=4 assumed
 
      xhh_arr_post_list_vb, phi_star_arr_list, a_star_arr_list, b_star_arr_list, Fphi_arr_list, Fa_arr_list, Fb_arr_list, logL_arr_list, F_arr_list, flag_converged_list, xhh_arr_post_list_all_try_vb, phi_star_arr_list_all_try, a_star_arr_list_all_try, b_star_arr_list_all_try, Fphi_arr_list_all_try, Fa_arr_list_all_try, Fb_arr_list_all_try, logL_arr_list_all_try, F_arr_list_all_try, flag_converged_list_all_try  = vbDCMM_global_func(
                                          xh_arr_post, p_xh_post, mu_arr_post,
                                          K_min=K_min, K_max=K_max,  
                                          uphi_prior=uphi_prior, 
                                          ua_prior=ua_prior, uad_prior=uad_prior,
                                          ub_prior=ub_prior, ubd_prior=ubd_prior,
                                          ua_prior_ic=ua_prior_ic, uad_prior_ic=uad_prior_ic,
                                          ub_prior_ic=ub_prior_ic, ubd_prior_ic=ubd_prior_ic,
                                          e_cutoff=1e-3, max_iter_vb=max_iter_vb, N_gl=N_gl,
                                          answer_ub_from_hmm=answer_ub_from_hmm)

      # Choose best model.
      _method = self.spinbox_F_methods.get()
      _Eend_arr = zeros( K_max - K_min + 1 ) # E = F or G
      _i = -1

      self.Fend_arr = zeros( self.K_max.get() - self.K_min.get() + 1 )
      self.Gend_arr = zeros( self.K_max.get() - self.K_min.get() + 1 )
      for k in range(K_min, K_max+1):
        _i += 1
        K_observed_set = set( xhh_arr_post_list_vb[_i] )
        Ko = len( K_observed_set )

        _F_arr = F_arr_list[_i]

        self.Fend_arr[_i] = _F_arr[-1] + log( gamma_func(k+1) ) # Model degeneracy.
        self.Gend_arr[_i] = _F_arr[-1] + log( gamma_func(Ko+1)*comb( k, Ko ) ) # Model degeneracy.
        if _method == 'Based on G':
          _Eend_arr[_i] = _F_arr[-1] + log( gamma_func(Ko+1)*comb( k, Ko ) ) # Model degeneracy.
        else: #'Based on F'
          _Eend_arr[_i] = _F_arr[-1] + log( gamma_func(k+1) ) # Model degeneracy.

      idx_max = argmax(_Eend_arr)
      _K = idx_max + K_min
      self.K_esti.set( _K )
      
      self.F_arr = _Eend_arr
      self.xhh_arr_VBDCMM = xhh_arr_post_list_vb[idx_max]
      self.ax_VBDCMM.hold(False)
      self.ax_VBDCMM.plot(self.time_xhh, self.xhh_arr_VBDCMM+1, 'k', linewidth=2)
      self.ax_VBDCMM.set_ylim([0.5, idx_max+1 + 0.5] )
      self.ax_VBDCMM.set_yticks( arange(1, idx_max+1+1, 1) )
      self.ax_VBDCMM.grid(True)
      self.ax_VBDCMM.hold(False)
      self.ax_VBDCMM.set_xlabel('Time (dt=1)')
      self.ax_VBDCMM.set_ylabel('X')
 

      self.fig_VBDCMM.canvas.draw()

      self.ax_F.hold(False)
      self.ax_F.plot(self.F_idx_arr, self.F_arr, '--rs')
      self.ax_F.grid(True)
      self.ax_F.hold(False)

      _xticks = arange( K_min, K_max+1, 1 )
      self.ax_F.set_xticks(_xticks )
      self.ax_F.set_xlim([K_min-0.5, K_max+0.5] )

      _yy = self.ax_F.get_ylim()
      _dyy = diff( _yy )
      _yy_min = _yy[0] - _dyy/7
      _yy_max = _yy[1] + _dyy/7
      self.ax_F.set_ylim([_yy_min, _yy_max] )
      self.ax_F.set_yticks(_yy)
      self.ax_F.set_xlabel('K')
      self.ax_F.set_ylabel('Evidence')

      self.fig_F.canvas.draw()

      # Update variables
      a_star_arr = a_star_arr_list[idx_max]
      a_star_arr_norm = normalize_transition_matrix( a_star_arr )
      str_exec = ''
      for i in range(1,self.K_esti_max+1):
        for j in range(1,self.K_esti_max+1):
          if i <= _K and j <= _K:
            str_exec += '\n'+'self.A' + str(i) + str(j) + '_esti.set(' + str( a_star_arr_norm[i-1, j-1] ) +')'
            str_exec += '\n'+'self.entry_A' + str(i) + str(j) + '_esti.config(relief=Tk.SOLID)'
          else:
            str_exec += '\n'+'self.A' + str(i) + str(j) + '_esti.set(' + str( 0 ) +')'
            str_exec += '\n'+'self.entry_A' + str(i) + str(j) + '_esti.config(relief=Tk.RIDGE)'
  
      exec(str_exec)

      b_star_arr = b_star_arr_list[idx_max]
      b_star_arr_norm = normalize_list_of_transition_matrix( b_star_arr )

      str_exec = ''
      for k in range(1, self.K_esti_max+1):
        for i in range(1,self.L_esti_max+1):
          for j in range(1,self.L_esti_max+1):
            # "exec" doesnt work inside of loop.
            if k <= _K and i <= _L and j <= _L:
              str_exec += '\n'+'self.B' + str(k) + str(i) + str(j) + '_esti.set(' + str( b_star_arr_norm[k-1,i-1,j-1] ) + ')'
              str_exec += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '_esti.config(relief=Tk.SOLID)'
            else:
              str_exec += '\n'+'self.B' + str(k) + str(i) + str(j) + '_esti.set(' + str( 0 ) + ')'
              str_exec += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '_esti.config(relief=Tk.RIDGE)'
      exec(str_exec)
    
 
  
    def on_K_change(self, *args):
      try:
        _K = self.K_esti.get()
        _L = self.L_esti.get()
        print( 'K: ' + str( _K ) )
        str_exec = ''
        for i in range(1, self.K_esti_max+1):
          for j in range(1, self.K_esti_max+1):
            #    if i >= _K+1 or j >= _K+1:
            #  str_exec += '\n'+'self.A' + str(i) + str(j) + '.set(0.0)'
            if i <= _K and j <= _K:
              str_exec += '\n'+'self.entry_A' + str(i) + str(j) + '.config(relief=Tk.SOLID)'
            else:
              str_exec += '\n'+'self.entry_A' + str(i) + str(j) + '.config(relief=Tk.RIDGE)'
  
        exec(str_exec)
        for k in range(1, self.K_esti_max+1):
          for i in range(1, self.L_esti_max+1):
            for j in range(1, self.L_esti_max+1):
              if k <= _K and i <= _L and j <= _L:
                str_exec += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '.config(relief=Tk.SOLID)'
              else:
                str_exec += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '.config(relief=Tk.RIDGE)'
        exec(str_exec)
  
        # special order
        if _K == 1:
          self.A11.set(1.0)
  
      except:
        pass
  
    def on_L_change(self, *args):
      try:
        _L = self.L_esti.get()
        _K = self.K_max.get()
        #print( 'N: ' + str( _L ) )
    
        str_exec = ''
        for i in range(1, self.L_esti_max+1):
          for j in range(1, self.L_esti_max+1):
            if i <= _L and j<= _L:
              str_exec += '\n'+'self.entry_B'  + str(i) + str(j)+'_hmm.config(relief=Tk.SOLID)'
              if i == j:
                str_exec += '\n'+'self.B' + str(i) + str(j) + '_hmm.set(' + str( 1 - (_L-1) * 0.01 ) + ')'
              else:
                str_exec += '\n'+'self.B' + str(i) + str(j) + '_hmm.set(' + str( 0.01 ) + ')'


            else:
              str_exec += '\n'+'self.entry_B'  + str(i) + str(j) + '_hmm.config(relief=Tk.RIDGE)'
              str_exec += '\n'+'self.B' + str(i) + str(j) + '_hmm.set(' + str( 0.0 ) + ')'
  
            for k in range(1, self.K_esti_max+1):
              if k <= _K and i <= _L and j <= _L:
                str_exec += '\n'+'self.entry_B' + str(k) + str(i) + str(j)+'_esti.config(relief=Tk.SOLID)'
              else:
                str_exec += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '_esti.config(relief=Tk.RIDGE)'
        exec(str_exec)
    
        str_exec = ''
        for i in range(1, self.L_esti_max+1):
          if i <= _L:
            str_exec += '\n'+'self.entry_mu' + str(i) + '_hmm.config(relief=Tk.SOLID)' 
            str_exec += '\n'+'self.entry_sig' + str(i) + '_hmm.config(relief=Tk.SOLID)' 
  
            str_exec += '\n'+'self.entry_mu' + str(i) + '_esti.config(relief=Tk.SOLID)' 
            str_exec += '\n'+'self.entry_sig' + str(i) + '_esti.config(relief=Tk.SOLID)' 
          else:
            str_exec += '\n'+'self.entry_mu' + str(i) + '_hmm.config(relief=Tk.RIDGE)' 
            str_exec += '\n'+'self.entry_sig' + str(i) + '_hmm.config(relief=Tk.RIDGE)' 
  
            str_exec += '\n'+'self.entry_mu' + str(i) + '_esti.config(relief=Tk.RIDGE)' 
            str_exec += '\n'+'self.entry_sig' + str(i) + '_esti.config(relief=Tk.RIDGE)' 
        exec(str_exec)

        _s = ''
        for i in range(1,_L+1):
          _s += '\n'+'self.mu'+str(i)+'_hmm.set(' + str( 0.1 + 0.8/(_L-1)* (i-1) ) + ')'
        exec(_s)
     
        _s = ''
        for i in range(1,_L+1):
          _s += '\n'+'self.sig'+str(i)+'_hmm.set('+ str( 0.15/_L ) + ')'
        exec(_s)
     

      except:
        #print('L exception!')
        pass

    #def on_B_hmm_change(self, *args):
    #  pass
  
    def on_mu_hmm_change(self, *args):
      pass
    def on_sig_hmm_change(self, *args):
      pass
  
    def on_T_xhh_change(self, *args):
      try:
        _T = self.T_xhh.get()
        self.time_xhh = arange(0, _T, 1)
        self.time_o = arange(0, _T+1, 1)
        print('T changed to %d.' % _T )
      except: # to catch error generated when the entry is empty.
        pass





  root = Tk.Tk()
  myapp = VBDCMM_analysis( root, o_arr )
  root.mainloop()

def generate_xhh_func( T_xhh, K, A ):
  """
  T_xhh: Length of trace
  K: The number of internal states
  A: Transition matrix for A

  return
  xhh_arr: Generated trace of internal states
  """
  # 1. Parameter & template
  T = T_xhh
  n_states_xhh = K # internal states
  
  # FRET value and noise
  p_xhh =deepcopy(A) # hidden state of hidden state
  
  xhh_arr = zeros( [T] )

  # xhh init. state
  r = rand() * n_states_xhh
  xhh_arr[0] = floor(r)

  # xhh states assign
  for i in range(1, T):
    r = rand()
    prev_n = int( xhh_arr[ i-1 ] )
    row_p_xhh = cumsum( p_xhh[prev_n,:] )

    # Assign next xhh state
    if r < row_p_xhh[0]:
      xhh_arr[i] = 0
    else:
      for i_xhh in range(1, n_states_xhh ):
        if r > row_p_xhh[ i_xhh - 1 ]:
          xhh_arr[i] = i_xhh
        else:
          break

  return xhh_arr
def generate_o_func( xhh_arr, K, L, B, mu_arr, sig_arr ):
  """
  T_xhh: Length of trace
  K: The number of internal states
  A: Transition matrix for A

  return
  xhh_arr: Generated trace of internal states
  """

  # 1. Parameter & template
  T = len( xhh_arr )
  
  n_states_xhh = K # internal states
  n_states_xh = L # Noise filtered observables
  
  # FRET value and noise
  vals_xh = mu_arr
  
  p_xh = deepcopy( B ) # hidden state
  
  xh_arr = zeros( [T+1] ) # xh_arr(1)-->xh_arr[1] affected by xhh_arr(1)
  o_arr = zeros( [T+1] )
  
  
  # xh, o init. state
  r = rand() * shape( p_xh )[1]
  xh_arr[0] = int( floor(r) )
  o_arr[0] = vals_xh[ int(xh_arr[0]) ] 
  
  row_p_xh = cumsum( p_xh[ int(xhh_arr[0]), int(xh_arr[0]), :] )
  r = rand()
  if r < row_p_xh[0]:
    xh_arr[1] = 0
    o_arr[1] = vals_xh[ 0 ] + randn(1)*sig_arr[0]
  else:
    for i_xh in range(1,  n_states_xh):
      if r > row_p_xh[ i_xh - 1 ]:
        xh_arr[1] = i_xh
        o_arr[1] = vals_xh[ i_xh ] + randn(1)*sig_arr[i_xh]
      else:
        break

  # xhh, xh states assign
  for i in range(1, T):
    r = rand()
    prev_n = xhh_arr[ i-1 ]
    prev_n_xh = xh_arr[ i-1 ]
    
    # Assign next xh state
    r = rand()
    row_p_xh = cumsum( p_xh[ int(xhh_arr[i]), int(xh_arr[i]), : ] )
    
    if r < row_p_xh[0]:
      xh_arr[i+1] = 0
      o_arr[i+1] = vals_xh[ 0 ] + randn(1)*sig_arr[0]
    else:
      for i_xh in range(1, n_states_xh):
        if r > row_p_xh[ i_xh - 1 ]:
          xh_arr[i+1] = i_xh
          o_arr[i+1] = vals_xh[ i_xh ] + randn(1)*sig_arr[i_xh]
        else:
          break

  return xh_arr, o_arr
