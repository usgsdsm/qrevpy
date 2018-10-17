import os
from PyQt5 import QtWidgets
import scipy.io as sio
# from Classes.stickysettings import StickySettings as SSet
import Panels.dsm_gui as dsm_gui
import sys
from Panels.selectFile import OpenMeasurementDialog, SaveMeasurementDialog
from Classes.Measurement import Measurement
# from Classes.ComputeExtrap import ComputeExtrap
# from Classes.QComp import QComp
from Classes.Python2Matlab import Python2Matlab
# from Classes.CommonDataComp import CommonDataComp


class TestDialog(QtWidgets.QMainWindow, dsm_gui.Ui_MainWindow):


    def __init__(self, parent=None):
        super(TestDialog, self).__init__(parent)
        self.settingsFile = 'TestSettings'
        # self.settings = SSet(self.settingsFile)
        self.setupUi(self)
        self.pushButton_4.clicked.connect(self.selectMeasurement)
        self.pushButton_3.clicked.connect(self.saveMeasurement)
    def selectMeasurement(self):
        self.select = OpenMeasurementDialog(self)
        self.select.exec_()
        if self.select.type == 'SonTek':
            # Show folder name in GUI header

            # Create measurement object
            self.meas = Measurement(in_file=self.select.fullName, source='SonTek', proc_type='QRev')

            print('SonTek')
        elif self.select.type == 'TRDI':
            # Show mmt filename in GUI header

            # Create measurement object
            self.meas = Measurement(in_file=self.select.fullName[0], source='TRDI', proc_type='QRev', checked=self.select.checked)
            self.meas.change_extrapolation(top='Constant', bot='NoSlip', exp=0.1)
            print('TRDI')
        elif self.select.type == 'QRev':
            self.meas = Measurement(in_file=self.select.fullName, source='QRev')
        else:
            print('Cancel')

    def saveMeasurement(self):
        # Create default file name
        save_file = SaveMeasurementDialog()
        Python2Matlab.save_matlab_file(self.meas, save_file.full_Name)




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