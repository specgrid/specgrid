#implements the Graphical user interface for

from qt_gui.mainwindow import Ui_MainWindow
from qt_gui.load_grid import Ui_Form
from mplwidget import MplWidget
from PyQt4 import QtCore, QtGui
import numpy as np

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure




class specgrid_mainwindow(QtGui.QMainWindow):
    def __init__(self, parent=None):
        self.setup_ui(parent)
        self.spec_plot = self.ax.plot([])[0]
        self.specgrids = []
        self.param_sliders = []
        self.autoscale = True
        self.off_grid_text = self.ax.text(0.5, 0.5,'OFF GRID',
                                horizontalalignment='center',
                                verticalalignment='center',
                                transform = self.ax.transAxes,
                                color='red',
                                fontsize=40)
    def setup_ui(self, parent):
        QtGui.QMainWindow.__init__(self, parent)
        self.main_frame = QtGui.QWidget(self)
        self.setCentralWidget(self.main_frame)
        self.resize(800, 600)
        
        self.vertical_layout = QtGui.QVBoxLayout()
        
        self.mplwidget = MplWidget(self.main_frame)
        self.ax = self.mplwidget.canvas.ax
        #self.vertical_layout.addWidget(self.mplwidget)
        self.autoscale_checkbox = QtGui.QCheckBox('Autoscale view')
        self.autoscale_checkbox.setChecked(True)
        self.autoscale_checkbox.stateChanged.connect(self.handle_autoscale_checkbox)
        self.vertical_layout.addWidget(self.mplwidget)
        self.vertical_layout.addWidget(self.autoscale_checkbox)
        
        self.main_frame.setLayout(self.vertical_layout)
        self.setWindowTitle('SpecGrid Mainwindow')

    def redraw(self):
        self.mplwidget.canvas.draw()
    
    def add_specgrid(self, specgrid):
        self.specgrids.append(specgrid)
        self.current_specgrid = specgrid
        self.current_grid_params = specgrid.param_mins.copy()
        new_flux = self.get_current_flux()
        
        self.test_slider = ParameterSlider(0, self)
        self.test_slider2 = ParameterSlider(1, self)
        
        for i in xrange(len(self.current_grid_params)):
            param_slider = ParameterSlider(i, self)
            self.param_sliders.append(param_slider)
            self.vertical_layout.addWidget(param_slider)
        self.plot_spec(new_flux, self.current_specgrid.wave, autoscale=True)
        
    def get_current_flux(self):
        flux = self.current_specgrid.interpolate(
                **dict(zip(self.current_specgrid.param_names, self.current_grid_params)))
        return flux
    
    def update_plot(self):
        flux = self.get_current_flux()
        self.plot_spec(flux)
    
    def plot_spec(self, flux, wave=None, autoscale=None):
        if all(np.isnan(flux)):
            self.off_grid_text.set_visible(True)
        else:
            self.off_grid_text.set_visible(False)
        
        if autoscale == None:
            autoscale = self.autoscale
        
        if wave!=None:
            self.spec_plot.set_xdata(wave)
        
        self.spec_plot.set_ydata(flux)
        
        if autoscale:
            self.ax.relim()
            self.ax.autoscale()
        
        self.redraw()
    
    def handle_autoscale_checkbox(self, checkbox_value):
        self.autoscale = bool(checkbox_value)

class ParameterSlider(QtGui.QWidget):
    def __init__(self, param_idx, main_window, parent=None):
        QtGui.QWidget.__init__(self, parent)
        
        self.param_idx = param_idx
        self.main_window = main_window
        self.current_specgrid = main_window.current_specgrid
        self.param_min = self.current_specgrid.param_mins[param_idx]
        self.param_max = self.current_specgrid.param_maxs[param_idx]
        self.param_name = self.current_specgrid.param_names[param_idx]
        self.current_params = self.main_window.current_grid_params
        
        
        self.horizontal_layout = QtGui.QHBoxLayout()
        
        self.param_label = QtGui.QLabel()
        self.param_label.setText(self.param_name)
        
        self.param_text_box = QtGui.QLineEdit()
        self.param_text_box.setText(str(self.current_params[param_idx]))
        self.param_text_box.editingFinished.connect(self.handle_text_box)
        self.param_linspace = np.linspace(self.param_min, self.param_max, 100)
        
        self.param_slider = QtGui.QSlider(1)
        self.param_slider.sliderMoved.connect(self.handle_slider)
        
        self.horizontal_layout.addWidget(self.param_label)
        self.horizontal_layout.addWidget(self.param_text_box)
        self.horizontal_layout.addWidget(self.param_slider)
        self.setLayout(self.horizontal_layout)
        
    def handle_slider(self, slider_value):
        new_param_value = self.param_linspace[slider_value]
        self.param_text_box.setText(str(new_param_value))
        self.current_params[self.param_idx] = new_param_value
        self.main_window.update_plot()
 
    def handle_text_box(self):
        try:
            new_param_value = np.float(self.param_text_box.text())
        except ValueError:
            self.param_text_box.setText(str(self.current_params[self.param_idx]))
            return
        
        new_param_value = min((max((new_param_value, self.param_min)), self.param_max))
        self.param_text_box.setText(str(new_param_value))
        self.param_slider.setValue(self.param_linspace.searchsorted(new_param_value))
        self.current_params[self.param_idx] = new_param_value
        self.main_window.update_plot()
 
        
class specgrid_load_grid(QtGui.QWidget):
    def __init__(self, parent=None):
        pass
        

if __name__ == '__main__':
    pass
