# 2015 Nov 7.
# Wonseok Hwang. 
# Load experimental data (1-D column)

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.pylab import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg 
from VBDCMM_gui_functions import VBDCMM_gui_analaysis_func

import matplotlib
#matplotlib.rcParams.update({'figure.autolayout': True})

import tkinter as Tk
from tkinter.filedialog import askopenfilenames

# 1. class
class VBDCMM_gui_load:
  def __init__( self, parent ):
    self.myParent = parent

    # Grids coordinates -----------------------------
    grid_button_load_data = (0, 2, 1, 2)
    grid_canvas_o = (0, 5, 1, 4*4+6)
    grid_button_analysis = (0, 30, 1, 2)

    # Variables ---------------------------
    self.load_filename = Tk.StringVar(master=parent, value='batch1')
#    self.label_load_filename = Tk.Label( master=parent, width=9, text='Filename:')
#    self.label_load_filename.grid( row=grid_label_load_filename[0], column=grid_label_load_filename[1], 
#        rowspan=grid_entry_load_filename[2],columnspan=grid_entry_load_filename[3])
#
#    self.entry_load_filename = Tk.Entry(master=parent, width=9, textvariable=self.load_filename )
#    self.entry_load_filename.grid( row=grid_entry_load_filename[0], column=grid_entry_load_filename[1],
#        rowspan=grid_entry_load_filename[2],columnspan=grid_entry_load_filename[3])
#
    # --- Button ---
    self.button_load_data = Tk.Button(text='Load data', width=8, command=self.on_click_button_load_data, height=4)
    self.button_load_data.grid( row=grid_button_load_data[0], column=grid_button_load_data[1],
                                        rowspan=grid_button_load_data[2], columnspan=grid_button_load_data[3],
                                        sticky='EW')
    self.button_analysis = Tk.Button(text='Analyze traces', width=10, command=self.on_click_button_analysis, height=4)
    self.button_analysis.grid( row=grid_button_analysis[0], column=grid_button_analysis[1],
                                        rowspan=grid_button_analysis[2], columnspan=grid_button_analysis[3],
                                        sticky='EW')


    # -------- Canvas --------
    _figsize_x = 7
    _figsize_y = 2
    _padx_1 = 0
    _pady_1 = 0
    _rect = [0.10, 0.3, 0.8, 0.6]
    self.o_arr = 0.5*ones(2200) + 0.1*randn( 2200 )
    self.time_o = arange(0, len(self.o_arr) )

    fign_o = 2

    self.fig_o = figure(num=fign_o, figsize=(_figsize_x, _figsize_y), facecolor='w')
    self.ax_o = self.fig_o.add_axes(_rect)
    self.ax_o.plot(self.time_o, self.o_arr, 'b')
    self.ax_o.set_ylim([-0.2, 1.2])
    self.ax_o.set_yticks([0, 0.5, 1.0])
    self.ax_o.grid(True)

    self.ax_o.hold(False)
 
    self.ax_o.set_xlabel('Time (dt=1)')
    self.ax_o.set_ylabel('O')

    self.canvas_o = FigureCanvasTkAgg(self.fig_o, master=parent)
    self.canvas_o.show()
    self.tk_canvas_o = self.canvas_o.get_tk_widget()
    self.tk_canvas_o.grid( row=grid_canvas_o[0], column=grid_canvas_o[1], 
                           rowspan=grid_canvas_o[2], columnspan=grid_canvas_o[3], 
                           padx=_padx_1, pady=_pady_1 )
    cid2 = self.canvas_o.mpl_connect('close_event', self.on_destroy)

  def on_destroy(self, event):
    """ 
    Without this..., loop after close the figure
    """
    print('Destroy!!') 
    self.myParent.quit()

  def on_click_button_load_data(self):
    #filename_0 = self.load_filename.get()
    filenames = askopenfilenames(initialdir='.') # it returns tuples of filename
    filename = filenames[0] 
    print('Load '+filename)
    self.o_arr = loadtxt( filename )
    _T = len(self.o_arr)
    self.time_o = arange(0, _T, 1)
    self.ax_o.plot(self.time_o, self.o_arr, 'b')
    _ymin = min( self.o_arr )
    _ymax = max( self.o_arr )

    self.ax_o.grid(True)
    self.ax_o.hold(False)
 
    self.ax_o.set_xlabel('Time (dt=1)')
    self.ax_o.set_ylabel('O')

    self.fig_o.canvas.draw()


  def on_click_button_analysis(self):
    print('VB-DCMM analysis')
    VBDCMM_gui_analaysis_func(self.o_arr)

    # 2. start
root = Tk.Tk()
myapp = VBDCMM_gui_load( root )

root.mainloop()
#Tk.mainloop()
# If you put root.destroy() here, it will cause an error if
# the window is closed with the window manager.
