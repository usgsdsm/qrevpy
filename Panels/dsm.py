import os
from PyQt5 import QtWidgets
import scipy.io as sio
from Classes.stickysettings import StickySettings as SSet
import dsm_gui
import sys
from Panels.selectFile import OpenMeasurementDialog
from Classes.Measurement import Measurement


class TestDialog(QtWidgets.QMainWindow, dsm_gui.Ui_MainWindow):


    def __init__(self, parent=None):
        super(TestDialog, self).__init__(parent)
        self.settingsFile = 'TestSettings'
        # self.settings = SSet(self.settingsFile)
        self.setupUi(self)
        self.pushButton_4.clicked.connect(self.selectMeasurement)

    def selectMeasurement(self):
        self.select = OpenMeasurementDialog(self)
        self.select.exec_()
        if self.select.type == 'SonTek':
            # Show folder name in GUI header

            # Create measurement object
            meas = Measurement(in_file=self.select.fullName, source='SonTek')

            print('SonTek')
        elif self.select.type == 'TRDI':
            # Show mmt filename in GUI header

            # Create measurement object
            meas = Measurement(in_file=self.select.fullName[0], source='TRDI', proc_type='QRev', checked=self.select.checked)
            print('TRDI')
        elif self.select.type == 'QRev':
            print('QRev')
        else:
            print('Cancel')




    def openFile(self):
        # Get the current folder setting. However, if the folder has not been previously defined create the Folder key
        # and set the value the current folder.
        try:
            folder = self.settings.get('Folder')
        except KeyError:
            self.settings.new('Folder',os.getcwd())
            folder = self.settings.get('Folder')
        # Allow the user to choose the file to open.
        fileName = QtWidgets.QFileDialog.getOpenFileName(self,self.tr('Open File'),folder,self.tr('Any File (*.mat)'))[0]
        # Update the folder setting
        pathName = os.path.split(fileName)[0]
        self.settings.set('Folder',pathName)
        # Read Matlab file
        mat_contents = sio.loadmat(fileName, struct_as_record=False, squeeze_me=True)
        # For QRev File
        meas_struct = mat_contents['meas_struct']
        QRev_version = mat_contents['version']
        print('Hello World')
        mat_struct = mat_contents['meas_struct']
app = QtWidgets.QApplication(sys.argv)
window = TestDialog()
window.show()
app.exec_()