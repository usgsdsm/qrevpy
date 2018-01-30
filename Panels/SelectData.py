from PyQt5.QtWidgets import (QWidget, QGridLayout, QPushButton, QCheckBox,
                             QLabel, QFileDialog)
from PyQt5.Qt import QRect, QDir
from functools import partial
import xmltodict

class SelectData(QWidget):
    """description of class"""

    def __init__(self):
        super(SelectData, self).init__()
        self.initUI()

    def initUI(self):
        self.layout = QGridLayout()
        self.setLayout(self.layout)
        self.dialog = QFileDialog()

        lbl = QLabel('<b style="font-size: 12pt">Select Options</b>')
        lbl.resize(250, 30)
        self.layout.addWidget(lbl,1,1,2,2)
        self.sontek = QPushButton("Sontek")
        self.sontek.clicked.connect(partial(self.uploadFile, uploadID = 1))
        self.trdi = QPushButton("TRDI")
        self.trdi.clicked.connect(partial(self.uploadFile, uploadID = 2))
        self.checkedOnly = QCheckBox()
        self.qrev = QPushButton("QRev")
        self.qrev.clicked.connect(partial(self.uploadFile, uploadID = 3))
        

        self.layout.addWidget(QLabel(''), 0,0)
        self.layout.addWidget(self.sontek, 2, 2)
        self.layout.addWidget(self.checkedOnly, 3, 1)
        self.layout.addWidget(self.trdi, 3, 2)
        self.layout.addWidget(self.qrev, 4, 2)
        self.layout.addWidget(QLabel(''), 5,3)

        self.layout.setColumnStretch(0, 10)
        self.layout.setColumnStretch(1, 1)
        self.layout.setColumnStretch(2, 1)
        self.layout.setColumnStretch(3, 10)
        self.layout.setRowStretch(0, 3)
        self.layout.setRowStretch(1, 1)
        self.layout.setRowStretch(2, 1)
        self.layout.setRowStretch(3, 1)
        self.layout.setRowStretch(4, 1)
        self.layout.setRowStretch(5, 5)
       
        

        self.setGeometry(QRect(0,0, 300, 300))

    def uploadFile(self, uploadID):

        if uploadID == 1:
                self.dialog.setNameFilter("*.mat")
        elif uploadID == 2:
                self.dialog.setNameFilter("*.mmt")

        self.dialog.show()

        
       
        
        


            

        

