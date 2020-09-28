
import sys
import os 
from PyQt5.QtWidgets import (QApplication,QMainWindow,QDialog,QFileDialog)

from qt_DWBA import MyWindow

def run_app():
        """
        launcher of Qt in Spyder
        to avoid error
        """
        if not QApplication.instance():
            app = QApplication(sys.argv)
        else:
            app = QApplication.instance()
        
        myWindow = MyWindow() 
        myWindow.show()
        sys.exit(app.exec_() )
        return myWindow
m = run_app()