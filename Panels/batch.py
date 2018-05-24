# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'batch.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_mainWindow(object):
    def setupUi(self, mainWindow):
        mainWindow.setObjectName("mainWindow")
        mainWindow.resize(266, 337)
        self.centralwidget = QtWidgets.QWidget(mainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.gb_create_list = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_create_list.setObjectName("gb_create_list")
        self.gridLayout = QtWidgets.QGridLayout(self.gb_create_list)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSpacing(12)
        self.verticalLayout.setObjectName("verticalLayout")
        self.rb_trdi = QtWidgets.QRadioButton(self.gb_create_list)
        self.rb_trdi.setChecked(True)
        self.rb_trdi.setObjectName("rb_trdi")
        self.verticalLayout.addWidget(self.rb_trdi)
        self.rb_sontek = QtWidgets.QRadioButton(self.gb_create_list)
        self.rb_sontek.setObjectName("rb_sontek")
        self.verticalLayout.addWidget(self.rb_sontek)
        self.rb_qrev = QtWidgets.QRadioButton(self.gb_create_list)
        self.rb_qrev.setObjectName("rb_qrev")
        self.verticalLayout.addWidget(self.rb_qrev)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 2, 1)
        self.pb_find = QtWidgets.QPushButton(self.gb_create_list)
        self.pb_find.setObjectName("pb_find")
        self.gridLayout.addWidget(self.pb_find, 0, 1, 1, 1)
        self.pb_save = QtWidgets.QPushButton(self.gb_create_list)
        self.pb_save.setObjectName("pb_save")
        self.gridLayout.addWidget(self.pb_save, 1, 1, 1, 1)
        self.gridLayout_2.addWidget(self.gb_create_list, 0, 0, 1, 1)
        self.gb_process = QtWidgets.QGroupBox(self.centralwidget)
        self.gb_process.setObjectName("gb_process")
        self.pb_process = QtWidgets.QPushButton(self.gb_process)
        self.pb_process.setGeometry(QtCore.QRect(20, 40, 75, 23))
        self.pb_process.setObjectName("pb_process")
        self.txt_status = QtWidgets.QLabel(self.gb_process)
        self.txt_status.setGeometry(QtCore.QRect(20, 90, 211, 41))
        self.txt_status.setText("")
        self.txt_status.setObjectName("txt_status")
        self.sb_skip = QtWidgets.QSpinBox(self.gb_process)
        self.sb_skip.setGeometry(QtCore.QRect(180, 40, 42, 22))
        self.sb_skip.setObjectName("sb_skip")
        self.label = QtWidgets.QLabel(self.gb_process)
        self.label.setGeometry(QtCore.QRect(150, 40, 21, 16))
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.gb_process, 1, 0, 1, 1)
        mainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(mainWindow)
        self.statusbar.setObjectName("statusbar")
        mainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(mainWindow)
        QtCore.QMetaObject.connectSlotsByName(mainWindow)

    def retranslateUi(self, mainWindow):
        _translate = QtCore.QCoreApplication.translate
        mainWindow.setWindowTitle(_translate("mainWindow", "QRevPy Batch Processing"))
        self.gb_create_list.setTitle(_translate("mainWindow", "Create List of Files"))
        self.rb_trdi.setText(_translate("mainWindow", "TRDI"))
        self.rb_sontek.setText(_translate("mainWindow", "SonTek"))
        self.rb_qrev.setText(_translate("mainWindow", "QRev"))
        self.pb_find.setText(_translate("mainWindow", "Find Files"))
        self.pb_save.setText(_translate("mainWindow", "Save List"))
        self.gb_process.setTitle(_translate("mainWindow", "Batch Process Files"))
        self.pb_process.setText(_translate("mainWindow", "Process Files"))
        self.label.setText(_translate("mainWindow", "Skip"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()
    ui = Ui_mainWindow()
    ui.setupUi(mainWindow)
    mainWindow.show()
    sys.exit(app.exec_())

