# 2015 Nov 7.
# Wonseok Hwang. 

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
from VBDCMM_gui_functions import generate_xhh_func
from VBDCMM_gui_functions import generate_o_func
from VBDCMM_gui_functions import VBDCMM_gui_analaysis_func

import matplotlib
#matplotlib.rcParams.update({'figure.autolayout': True})

import tkinter as Tk

# 1. class
class VBDCMM_gui_simul:
  def __init__( self, parent ):
    self.myParent = parent

    # --- Model parameters ------------------------------
    self.K_true = Tk.IntVar(master=parent)
    self.K_true.trace('w', self.on_K_change)
    self.K_true_max = 4
    self.L_true = Tk.IntVar(master=parent)
    self.L_true.trace('w', self.on_L_change )
    self.L_true_max = 4

    self.T_xhh = Tk.IntVar(master=parent)
    self.T_xhh.trace('w', self.on_T_xhh_change )
    self.T_o = Tk.IntVar(master=parent)
    self.T_xhh_default = 4400
    self.T_o_default = self.T_xhh_default+1

    self.time_xhh = arange(0, self.T_xhh_default, 1)
    self.xhh_arr = zeros( self.T_xhh_default )  # internal states index starts from 0. However, for plotting, it starts from 1
    self.time_o = arange(0, self.T_o_default, 1) 
    self.xh_arr = zeros( self.T_o_default )  # internal states index starts from 0. However, for plotting, it starts from 1
    self.o_arr = 0.5*ones(self.T_o_default ) + 0.1*randn( self.T_o_default )

    # Grids coordinates -----------------------------
    grid_label_K_true = (0, 0)  # row, column, rowspan, columnspan
    grid_entry_K_true = (0, 1)
    grid_label_T_xhh = (0, 2)
    grid_entry_T_xhh = (0, 3)
    grid_label_A = (0, 4)
    grid_entry_A = (0, 5)
    grid_button_generate_params = (1, 5+self.K_true_max, 2, 2)
    grid_button_generate_trace_xhh = (0+self.K_true_max, 4, 1, 2)
    grid_canvas_xhh = (0+self.K_true_max, 0, 1, self.L_true_max*4 + 6)

    grid_label_L_true = (0+self.K_true_max+1, 0)  # row, column, rowspan, columnspan
    grid_entry_L_true = (0+self.K_true_max+1, 1)
    grid_label_mu_true = (0+self.K_true_max+1, 2)
    grid_entry_mu1_true = (0+self.K_true_max+1, 3)
    grid_entry_mu2_true = (0+self.K_true_max+2, 3)
    grid_entry_mu3_true = (0+self.K_true_max+3, 3)
    grid_entry_mu4_true = (0+self.K_true_max+4, 3)
    grid_label_sig_true = (0+self.K_true_max+1, 4)
    grid_entry_sig1_true = (0+self.K_true_max+1, 5)
    grid_entry_sig2_true = (0+self.K_true_max+2, 5)
    grid_entry_sig3_true = (0+self.K_true_max+3, 5)
    grid_entry_sig4_true = (0+self.K_true_max+4, 5)

    grid_label_B1 = (0+self.K_true_max+1, 6)
    grid_entry_B1 = (0+self.K_true_max+1, 7)
    grid_label_B2 = (0+self.K_true_max+1, 7+self.L_true_max)
    grid_entry_B2 = (0+self.K_true_max+1, 7+self.L_true_max+1)
    grid_label_B3 = (0+self.K_true_max+1, 7+2*self.L_true_max+1)
    grid_entry_B3 = (0+self.K_true_max+1, 7+2*self.L_true_max+2)
    grid_label_B4 = (0+self.K_true_max+1, 7+3*self.L_true_max+2)
    grid_entry_B4 = (0+self.K_true_max+1, 7+3*self.L_true_max+3)
 
    grid_button_generate_trace_o = (0+self.K_true_max+1+self.L_true_max, 4, 1, 2)
    grid_canvas_o = (0+self.K_true_max+1+self.L_true_max, 0, 1, self.L_true_max*4+6)

    grid_label_save_filename = (0+1*self.K_true_max+1+self.L_true_max, 17, 1, 2)
    grid_entry_save_filename = (0+1*self.K_true_max+1+self.L_true_max, 19, 1, 2)
    grid_button_save_traces = (0+1*self.K_true_max+1+self.L_true_max, 22, 1, 2)
    grid_button_analysis = (0+1*self.K_true_max+1+self.L_true_max, 24, 1, 2)

    # Variables ---------------------------
    _s = ''
    for i in range(1,self.K_true_max+1):
      for j in range(1,self.K_true_max+1):
        _s += '\n'+'self.A' + str(i) + str(j) + ' = Tk.DoubleVar(master=parent)'
        _s += '\n'+'self.A' + str(i) + str(j) + '.trace("w", self.on_A_change)'
    exec(_s)

    _s = ''
    for k in range(1, self.K_true_max+1):
      for i in range(1,self.L_true_max+1):
        for j in range(1,self.L_true_max+1):
          _s += '\n'+'self.B' +str(k)  + str(i) + str(j) + ' = Tk.DoubleVar(master=parent)'
          _s += '\n'+'self.B' +str(k)  + str(i) + str(j) + '.trace("w", self.on_B_change)'
    exec(_s)
 
    _s = ''
    for i in range(1,self.L_true_max+1):
      _s += '\n'+'self.mu'+str(i)+'_true = Tk.DoubleVar(master=parent)'
      _s += '\n'+'self.mu'+str(i)+'_true.trace("w", self.on_mu_change)'
      _s += '\n'+'self.sig'+str(i)+'_true = Tk.DoubleVar(master=parent)'
      _s += '\n'+'self.sig'+str(i)+'_true.trace("w", self.on_sig_change)'
    exec(_s)
       
      
    # --- Set initial values -----------------------------
    self.K_default = 2
    self.K_true.set(self.K_default)
    self.T_xhh.set(self.T_xhh_default)
    self.T_o.set(self.T_o_default)
    self.L_default = 2
    self.L_true.set(self.L_default)

    _s = ''
    for i in range(1,self.K_default+1):
      for j in range(1,self.K_default+1):
        if i != j:
          _s += '\n'+'self.A' + str(i) + str(j) + '.set(' + str( 0.001/(self.K_default-1)) + ')'
        if i == j:
          _s += '\n'+'self.A' + str(i) + str(j) + '.set(' + str( 1 - 0.001) + ')'
    exec(_s)
 
    _s = ''
    for k in range(1, self.K_default+1):
      for i in range(1,self.L_default+1):
        for j in range(1,self.L_default+1):
          if i != j:
            _s += '\n'+'self.B' + str(k) + str(i) + str(j) + '.set(' + str( (1/4**(k-1) * 0.1)/(self.L_default-1)) + ')'
          if i == j:
            _s += '\n'+'self.B' + str(k) + str(i) + str(j) + '.set(' + str( 1 - 1/4**(k-1) * 0.1) + ')'
    exec(_s)

    _s = ''
    for i in range(1,self.L_default+1):
      _s += '\n'+'self.mu'+str(i)+'_true.set(' + str( 0.1 + 0.8/(self.L_default-1)* (i-1) ) + ')'
    exec(_s)
 
    _s = ''
    for i in range(1,self.L_default+1):
      _s += '\n'+'self.sig'+str(i)+'_true.set('+ str( 0.15/self.L_default ) + ')'
    exec(_s)
 
    # --- Set trace function --- 
    #self.K_true.trace('w', on_change_K_true)


    # --- Labels ---
    self.label_K_true = Tk.Label( master=parent, width=2, text='K:' )
    self.label_A = Tk.Label( master=parent, width=2, text='A:' )
    self.label_T_xhh = Tk.Label( master=parent, width=2, text='T:' )
    self.label_L_true = Tk.Label( master=parent, width=2, text='N:' )
    self.label_mu_true = Tk.Label( master=parent, width=2, text='mu:')
    self.label_sig_true = Tk.Label( master=parent, width=2, text='sig:')
    self.label_B1 = Tk.Label( master=parent, width=2, text='B1:')
    self.label_B2 = Tk.Label( master=parent, width=2, text='B2:')
    self.label_B3 = Tk.Label( master=parent, width=2, text='B3:')
    self.label_B4 = Tk.Label( master=parent, width=2, text='B4:')

    # --- Labels grid ---
    self.label_K_true.grid( row=grid_label_K_true[0], column=grid_label_K_true[1] )
    self.label_A.grid( row=grid_label_A[0], column=grid_label_A[1] )
    self.label_T_xhh.grid( row=grid_label_T_xhh[0], column=grid_label_T_xhh[1] )
    self.label_L_true.grid( row=grid_label_L_true[0], column=grid_label_L_true[1] )
    self.label_mu_true.grid( row=grid_label_mu_true[0], column=grid_label_mu_true[1] )
    self.label_sig_true.grid( row=grid_label_sig_true[0], column=grid_label_sig_true[1] )
    self.label_B1.grid( row=grid_label_B1[0], column=grid_label_B1[1] )
    self.label_B2.grid( row=grid_label_B2[0], column=grid_label_B2[1] )
    self.label_B3.grid( row=grid_label_B3[0], column=grid_label_B3[1] )
    self.label_B4.grid( row=grid_label_B4[0], column=grid_label_B4[1] )

    self.save_filename = Tk.StringVar(master=parent, value='batch1')
    self.label_save_filename = Tk.Label( master=parent, width=9, text='Filename:')
    self.label_save_filename.grid( row=grid_label_save_filename[0], column=grid_label_save_filename[1], 
        rowspan=grid_entry_save_filename[2],columnspan=grid_entry_save_filename[3])

    self.entry_save_filename = Tk.Entry(master=parent, width=9, textvariable=self.save_filename )
    self.entry_save_filename.grid( row=grid_entry_save_filename[0], column=grid_entry_save_filename[1],
        rowspan=grid_entry_save_filename[2],columnspan=grid_entry_save_filename[3])


    # --- Button ---
    self.button_generate_params = Tk.Button(text='Gen. Params', width=8, command=self.on_click_button_generate_params, height=4)
    self.button_generate_params.grid( row=grid_button_generate_params[0], column=grid_button_generate_params[1],
                                        rowspan=grid_button_generate_params[2], columnspan=grid_button_generate_params[3],
                                        sticky='EW')
    self.button_save_traces = Tk.Button(text='Save. traces', width=10, command=self.on_click_button_save_traces, height=4)
    self.button_save_traces.grid( row=grid_button_save_traces[0], column=grid_button_save_traces[1],
                                        rowspan=grid_button_save_traces[2], columnspan=grid_button_save_traces[3],
                                        sticky='EW')
    self.button_analysis = Tk.Button(text='Analyze traces', width=8, command=self.on_click_button_analysis, height=4)
    self.button_analysis.grid( row=grid_button_analysis[0], column=grid_button_analysis[1],
                                        rowspan=grid_button_analysis[2], columnspan=grid_button_analysis[3],
                                        sticky='EW')



    self.button_generate_trace_xhh = Tk.Button(text='Gen. X', width=8, command=self.on_click_button_generate_trace_xhh, relief=Tk.SUNKEN, height=4)
    self.button_generate_trace_xhh.grid( row=grid_button_generate_trace_xhh[0], column=grid_button_generate_trace_xhh[1],
                                        rowspan=grid_button_generate_trace_xhh[2], columnspan=grid_button_generate_trace_xhh[3],
                                        sticky='EW')

    self.button_generate_trace_o = Tk.Button(text='Gen. O', width=8, command=self.on_click_button_generate_trace_o, relief=Tk.SUNKEN, height=4)
    self.button_generate_trace_o.grid( row=grid_button_generate_trace_o[0], column=grid_button_generate_trace_o[1],
                                        rowspan=grid_button_generate_trace_o[2], columnspan=grid_button_generate_trace_o[3],
                                        sticky='EW')


    # --- Entries --- 
    #self.entry_K_true = Tk.Entry( relief=Tk.SOLID, master=parent, textvariable=self.K_true, width=3) # SOLID: style. Clearer boundary
    self.entry_K_true = Tk.Spinbox( relief=Tk.RIDGE, master=parent, textvariable=self.K_true, width=3, value=(1,2,3,4) ) # SOLID: style. Clearer boundary
    self.K_true.set(self.K_default)
    self.entry_T_xhh = Tk.Entry( relief=Tk.SOLID, master=parent, textvariable=self.T_xhh, width=5) # SOLID: style. Clearer boundary
    #self.entry_L_true = Tk.Entry( relief=Tk.SOLID, master=parent, textvariable=self.L_true, width=3) # SOLID: style. Clearer boundary
    self.entry_L_true = Tk.Spinbox( relief=Tk.RIDGE, master=parent, textvariable=self.L_true, width=3, value=(1,2,3,4) ) # SOLID: style. Clearer boundary
    self.L_true.set(self.L_default)
    _s = ''
    for i in range(1,self.K_true_max+1):
      for j in range(1,self.K_true_max+1):
        _s += '\n'+'self.entry_A' + str(i) + str(j) + ' = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.A'+str(i)+str(j)+', width=7)'
        if i <= self.K_default and j <= self.K_default:
          _s += '\n'+'self.entry_A' + str(i) + str(j) + '.config(relief=Tk.SOLID)'
    exec(_s)

    _s = ''
    for i in range(1,self.L_true_max+1):
      _s += '\n'+'self.entry_mu' + str(i) + '_true = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.mu'+str(i)+'_true, width=5)'
      _s += '\n'+'self.entry_sig' + str(i) + '_true = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.sig'+str(i)+'_true, width=5)'
      if i <= self.L_default:
        _s += '\n'+'self.entry_mu' + str(i) + '_true.config(relief=Tk.SOLID)' 
        _s += '\n'+'self.entry_sig' + str(i) + '_true.config(relief=Tk.SOLID)' 
    exec(_s)

    _s = ''
    for k in range(1, self.K_true_max+1): 
      for i in range(1,self.L_true_max+1):
        for j in range(1,self.L_true_max+1):
          _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + ' = Tk.Entry(relief=Tk.RIDGE, master=parent, textvariable=self.B'+str(k)+str(i)+str(j)+', width=5)'
          if k <= self.K_default and i <= self.L_default and j <=self.L_default:
            _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '.config(relief=Tk.SOLID)'
    exec(_s)
   

    # --- Entry grid --- 
    self.entry_K_true.grid(row=grid_entry_K_true[0], column=grid_entry_K_true[1], pady=0.1)
    self.entry_T_xhh.grid(row=grid_entry_T_xhh[0], column=grid_entry_T_xhh[1], pady=0.1)
    self.entry_L_true.grid(row=grid_entry_L_true[0], column=grid_entry_L_true[1], pady=0.1)

    _s = ''
    for i in range(1,self.K_true_max+1):
      for j in range(1,self.K_true_max+1):
        _row = (i-1) + grid_entry_A[0]
        _column = (j-1) + grid_entry_A[1]
        _s += '\n'+'self.entry_A' + str(i) + str(j) + '.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
    exec(_s)
 
    _s = ''
    for i in range(1,self.L_true_max+1):
      _row = eval( 'grid_entry_mu' + str(i) + '_true[0]')
      _column = eval( 'grid_entry_mu' + str(i) + '_true[1]')
      _s += '\n'+'self.entry_mu' + str(i) + '_true.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
    exec(_s)

    _s = ''
    for i in range(1,self.L_true_max+1):
      _row = eval( 'grid_entry_sig' + str(i) + '_true[0]')
      _column = eval( 'grid_entry_sig' + str(i) + '_true[1]')
      _s += '\n'+'self.entry_sig' + str(i) + '_true.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
    exec(_s)
 
    _s = ''
    for k in range(1, self.K_true_max+1):
      for i in range(1,self.L_true_max+1):
        for j in range(1,self.L_true_max+1):
          # "exec" doesnt work inside of loop.
          _row = (i-1) + eval('grid_entry_B'+str(k)+'[0]' )
          _column = (j-1) + eval('grid_entry_B'+str(k)+'[1]' )

          _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '.grid(row='+str(_row)+', column='+str(_column)+', pady=0.1)'    
    exec(_s)
   
  
    # -------- Canvas --------
    _figsize_x = 7
    _figsize_y = 2
    _padx_1 = 0
    _pady_1 = 0
    _rect = [0.10, 0.3, 0.8, 0.6]
    fign_xhh = 1
    fign_o = 2
    fign_xhh_VBDCMM = 3

    self.fig_xhh = figure(num=fign_xhh, figsize=(_figsize_x, _figsize_y), facecolor='w')
    self.ax_xhh = self.fig_xhh.add_axes(_rect)
    self.ax_xhh.plot(self.time_xhh, self.xhh_arr+1,'k', linewidth=2)
    self.ax_xhh.set_ylim([0.5, self.K_default + 0.5] )
    self.ax_xhh.set_yticks( arange(1, self.K_default+1, 1) )
    self.ax_xhh.grid(True)
 
    self.ax_xhh.hold(False)

    self.canvas_xhh = FigureCanvasTkAgg(self.fig_xhh, master=parent)
    self.canvas_xhh.show()
    self.tk_canvas_xhh = self.canvas_xhh.get_tk_widget()
    self.tk_canvas_xhh.grid( row=grid_canvas_xhh[0], column=grid_canvas_xhh[1], 
                             rowspan=grid_canvas_xhh[2], columnspan=grid_canvas_xhh[3], 
                             padx=_padx_1, pady=_pady_1 )
 
    xlabel('Time (dt=1)')
    ylabel('X')

    cid1 = self.canvas_xhh.mpl_connect('button_press_event', self.on_figure_click)
    cid2 = self.canvas_xhh.mpl_connect('close_event', self.on_destroy)

    self.fig_o = figure(num=fign_o, figsize=(_figsize_x, _figsize_y), facecolor='w')
    self.ax_o = self.fig_o.add_axes(_rect)
    self.ax_o.plot(self.time_o, self.o_arr,'b')
    self.ax_o.set_ylim([-0.2, 1.2])
    self.ax_o.set_yticks([0, 0.5, 1.0])
    self.ax_o.grid(True)

    self.ax_o.hold(False)

    self.canvas_o = FigureCanvasTkAgg(self.fig_o, master=parent)
    self.canvas_o.show()
    self.tk_canvas_o = self.canvas_o.get_tk_widget()
    self.tk_canvas_o.grid( row=grid_canvas_o[0], column=grid_canvas_o[1], 
                           rowspan=grid_canvas_o[2], columnspan=grid_canvas_o[3], 
                           padx=_padx_1, pady=_pady_1 )
 
    xlabel('Time (dt=1)')
    ylabel('O')

    cid1 = self.canvas_o.mpl_connect('button_press_event', self.on_figure_click)
    cid2 = self.canvas_o.mpl_connect('close_event', self.on_destroy)

  def on_destroy(self, event):
    """ 
    Without this..., loop after close the figure
    """
    print('Destroy!!') 
    self.myParent.quit()

  def on_figure_click(self, event):
    print('you pressed', event.button, event.x, event.y, event.xdata, event.ydata)
    #print('e1_number: ' + str(self.e1_number.get()) )

  def on_click_button_generate_params(self):
    print('Generated parameters')
    _K = self.K_true.get()
    _L = self.L_true.get()
    _T = self.T_xhh.get()

    _s = ''
    for i in range(1,_K+1):
      for j in range(1,_K+1):
        if i != j:
          _s += '\n'+'self.A' + str(i) + str(j) + '.set(' + str( 0.0005/(_K-1)) + ')'
        if i == j:
          _s += '\n'+'self.A' + str(i) + str(j) + '.set(' + str( 1 - 0.0005) + ')'
    exec(_s)
 
    _s = ''
    for k in range(1, _K+1):
      for i in range(1, _L+1):
        for j in range(1, _L+1):
          if i != j:
            _s += '\n'+'self.B' + str(k) + str(i) + str(j) + '.set(' + str( (1/4**(k-1) * 0.1)/(_L-1)) + ')'
          if i == j:
            _s += '\n'+'self.B' + str(k) + str(i) + str(j) + '.set(' + str( 1 - 1/4**(k-1) * 0.1) + ')'

    exec(_s)

    _s = ''
    for i in range(1,_L+1):
      _s += '\n'+'self.mu'+str(i)+'_true.set(' + str( 0.1 + 0.8/(_L-1)* (i-1) ) + ')'
    exec(_s)
 
    _s = ''
    for i in range(1,_L+1):
      _s += '\n'+'self.sig'+str(i)+'_true.set('+ str( 0.15/_L ) + ')'
    exec(_s)
 
  def on_click_button_save_traces(self):
    print('Save traces to %s.dat' % self.save_filename.get())
    filename_0 = self.save_filename.get()
    savetxt(filename_0 + '_X_true_trace.dat', self.xhh_arr, fmt='%d')
    savetxt(filename_0 + '_O_filt_true.dat', self.xh_arr, fmt='%.3f' )
    savetxt(filename_0 + '_O_true.dat', self.o_arr, fmt='%.3f' )

    # save A
    _K = self.K_true.get()
    _T = self.T_xhh.get()
    self.A = zeros( [_K, _K] ) 
    for i in range(1, _K+1):
      for j in range(1, _K+1):
        self.A[i-1, j-1] = eval( 'self.A' + str(i) + str(j) + '.get()' )
    savetxt(filename_0+'_A_true.dat', self.A, fmt='%.9f' )

    _L = self.L_true.get()
    self.B = zeros( [_K, _L, _L] ) 
    self.mu_arr = zeros( _L )
    self.sig_arr = zeros( _L )

    # save B
    for k in range(1, _K+1):
      for i in range(1, _L+1):
        for j in range(1, _L+1):
          self.B[k-1, i-1, j-1] = eval( 'self.B' + str(k) + str(i) + str(j) + '.get()' )

      savetxt(filename_0+'_B'+str(k)+'_true.dat', self.B[k-1,:,:], fmt='%.9f' )

    # save mu and sig
    self.mu_arr = zeros( _L )
    self.sig_arr = zeros( _L )

    for i in range(1, _L+1):
      self.mu_arr[i-1] = eval('self.mu'+str(i)+'_true.get()')
      self.sig_arr[i-1] = eval('self.sig'+str(i)+'_true.get()')

    savetxt(filename_0+'_mu_arr_true.dat', self.mu_arr, fmt='%.3f' )
    savetxt(filename_0+'_sig_arr_true.dat', self.sig_arr, fmt='%.3f' )

 
  def on_click_button_analysis(self):
    print('VB-DCMM analysis')
    VBDCMM_gui_analaysis_func(self.o_arr)

  def on_click_button_generate_trace_xhh(self):
    print('Generate trace_xhh')
    # Construct A matrix
    _K = self.K_true.get()
    _T = self.T_xhh.get()
    self.A = zeros( [_K, _K] ) 
    for i in range(1, _K+1):
      for j in range(1, _K+1):
        self.A[i-1, j-1] = eval( 'self.A' + str(i) + str(j) + '.get()' )

    self.xhh_arr = generate_xhh_func( _T, _K, self.A ) 
    self.ax_xhh.plot(self.time_xhh, self.xhh_arr+1, 'k', linewidth=2)
    self.ax_xhh.set_ylim([0.5, _K + 0.5] )
    self.ax_xhh.set_yticks( arange(1, _K+1, 1) )
    self.ax_xhh.grid(True)
    self.ax_xhh.set_xlabel('Time (dt=1)')
    self.ax_xhh.set_ylabel('X')
    self.fig_xhh.canvas.draw()

  def on_click_button_generate_trace_o(self):
    print('Generate trace_xhh')
    # Construct A matrix
    _K = self.K_true.get()
    _L = self.L_true.get()
    _T = self.T_xhh.get()
    self.A = zeros( [_K, _K] ) 
    self.B = zeros( [_K, _L, _L] ) 
    self.mu_arr = zeros( _L )
    self.sig_arr = zeros( _L )

    for k in range(1, _K+1):
      for i in range(1, _L+1):
        for j in range(1, _L+1):
          self.B[k-1, i-1, j-1] = eval( 'self.B' + str(k) + str(i) + str(j) + '.get()' )

    for i in range(1, _L+1):
      self.mu_arr[i-1] = eval('self.mu'+str(i)+'_true.get()')
      self.sig_arr[i-1] = eval('self.sig'+str(i)+'_true.get()')

    self.xh_arr, self.o_arr = generate_o_func( self.xhh_arr, _K, _L, self.B, self.mu_arr, self.sig_arr ) 
    self.ax_o.plot(self.time_o, self.o_arr)
    self.ax_o.set_ylim([-0.2, 1.2])
    self.ax_o.set_yticks([0, 0.5, 1.0])
    self.ax_o.grid(True)
    self.ax_o.set_xlabel('Time (dt=1)')
    self.ax_o.set_ylabel('O')
    self.fig_o.canvas.draw()

  def on_K_change(self, *args):
    try:
      _K = self.K_true.get()
      _L = self.L_true.get()
      #print( 'K: ' + str( _K ) )
      _s = ''
      for i in range(1, self.K_true_max+1):
        for j in range(1, self.K_true_max+1):
          #    if i >= _K+1 or j >= _K+1:
          #  _s += '\n'+'self.A' + str(i) + str(j) + '.set(0.0)'
          if i <= _K and j <= _K:
            _s += '\n'+'self.entry_A' + str(i) + str(j) + '.config(relief=Tk.SOLID)'
          else:
            _s += '\n'+'self.entry_A' + str(i) + str(j) + '.config(relief=Tk.RIDGE)'

      exec(_s)
      for k in range(1, self.K_true_max+1):
        for i in range(1, self.L_true_max+1):
          for j in range(1, self.L_true_max+1):
            if k <= _K and i <= _L and j <= _L:
              _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '.config(relief=Tk.SOLID)'
            else:
              _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '.config(relief=Tk.RIDGE)'
      exec(_s)

      # special order
      if _K == 1:
        self.A11.set(1.0)

    except:
      pass

  def on_L_change(self, *args):
    try:
      _K = self.K_true.get()
      _L = self.L_true.get()
      #print( 'N: ' + str( _L ) )

      _s = ''
      for k in range(1, self.K_true_max+1):
        for i in range(1, self.L_true_max+1):
          for j in range(1, self.L_true_max+1):
            if k <= _K and i <= _L and j <= _L:
              _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j)+'.config(relief=Tk.SOLID)'
            else:
              _s += '\n'+'self.entry_B' + str(k) + str(i) + str(j) + '.config(relief=Tk.RIDGE)'
      exec(_s)

      _s = ''
      for i in range(1, self.L_true_max+1):
        if i <= _L:
          _s += '\n'+'self.entry_mu' + str(i) + '_true.config(relief=Tk.SOLID)' 
          _s += '\n'+'self.entry_sig' + str(i) + '_true.config(relief=Tk.SOLID)' 
        else:
          _s += '\n'+'self.entry_mu' + str(i) + '_true.config(relief=Tk.RIDGE)' 
          _s += '\n'+'self.entry_sig' + str(i) + '_true.config(relief=Tk.RIDGE)' 
      exec(_s)
    except:
      pass
 


  def on_A_change(self, *args):
      pass

  def on_B_change(self, *args):
    pass

  def on_mu_change(self, *args):
    pass
  def on_sig_change(self, *args):
    pass

  def on_T_xhh_change(self, *args):
    try:
      _T = self.T_xhh.get()
      self.time_xhh = arange(0, _T, 1)
      self.time_o = arange(0, _T+1, 1)
      #print('T changed to %d.' % _T )
    except: # to catch error generated when the entry is empty.
      pass

    # 2. start
root = Tk.Tk()
myapp = VBDCMM_gui_simul( root )

root.mainloop()
