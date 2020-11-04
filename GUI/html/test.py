import sys
import os
from PyQt5.QtWidgets import (QApplication,QMainWindow,QDialog,QFileDialog,QWidget,
                             QComboBox,QLabel,QLineEdit,QCheckBox,
                             QHBoxLayout,QVBoxLayout,QGridLayout,
                             QGroupBox,QToolBox,QTabWidget,QPushButton,QTextBrowser,
                             QSpacerItem,QRadioButton)

from PyQt5 import uic
from PyQt5 import QtCore

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas,
        NavigationToolbar2QT as NavigationToolbar)

from numpy import loadtxt

from io import StringIO
from subprocess import (call,Popen)

import numpy as np


class text_Browser(QWidget,):
    """
    Change between Text Browser and Web Browser
    """
    def __init__(self,browser_choice='text',initial_page=''):
        """

        browser_choice =='text' for text_browser
                               but this can also show simple Html text
                       =='web'  for web browser

        set_page_txt puts text either plain string text
        or HTML format

        """
        from PyQt5.QtWebEngineWidgets import QWebEngineView
        from PyQt5.QtCore import QUrl

        super().__init__()
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.browser_choice = browser_choice
        self.page_txt = initial_page
        if self.browser_choice =='text':
            self.browser = QTextBrowser()
            self.browser.setAcceptRichText(True)
            self.browser.setOpenExternalLinks(True)
        elif self.browser_choice =='web':
            self.browser = QWebEngineView()
            pass
        #---initial view
        self.layout.addWidget(self.browser)
        self.set_page_txt(self.page_txt)

    def set_page_txt(self,page_txt,txt_type='text'):
        self.page_txt = page_txt
        if self.browser_choice =='text':
            if txt_type=='text':
                self.browser.setPlainText(page_txt)
            elif txt_type=='html':
                self.browser.setHtml(page_txt)
        elif self.browser_choice =='web':
            self.browser.setHtml(page_txt)
            
if __name__ == "__main__":

  import sys
  from PyQt5.Qt import *
  from PyQt5.QtWebEngineWidgets import *
  from PyQt5.QtWidgets import QApplication
 
  f=open('description_transfer.html','r')
  hhh = f.read()  
  f.close()
 
  app = QApplication(sys.argv)
 
  web = QWebEngineView()
 
  #web.load(QUrl("https://www.codeloop.org"))
  web.setHtml(hhh)
 
  web.show()
 
 
  sys.exit(app.exec_())
    




 