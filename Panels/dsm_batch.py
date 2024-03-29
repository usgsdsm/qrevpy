import os
import sys
import copy
from PyQt5 import QtWidgets
from PyQt5 import QtCore
import Panels.batch as batch
from Classes.Measurement import Measurement
import xlsxwriter
import numpy as np


class BatchDialog(QtWidgets.QMainWindow, batch.Ui_mainWindow):
    """Allows the user to create a file contain a list of all TRDI mmt, SonTek mat, or QRev mat files starting at
    the specified top level folder and recursively searching all folders below. The option is then presented to
    batch process these files using QRev.

    Attributes
    ----------
    files: list
        List of files with full path meeting requested criteria.
    top_folder: str
        User selected top folder to begin searching.
    """

    def __init__(self, parent=None):
        """Initializes settings and connections.

        Parameters
        ----------
        parent
            Identifies parent GUI.
        """

        super(BatchDialog, self).__init__(parent)
        self.setupUi(self)

        # Initialize variables
        self.files = []
        self.top_folder = ''
        self.pb_save.setEnabled(False)

        # Create connections for buttons
        self.pb_find.clicked.connect(self.find_files)
        self.pb_save.clicked.connect(self.save_list)
        self.pb_process.clicked.connect(self.process)

        self.summary = []
    def find_files(self):
        """Finds all files in the folder and subfolder matching the criteria and saves them to a list."""

        self.files = []
        # Disable save button until list is created
        self.pb_save.setEnabled(False)

        # Get top level folder from user
        self.top_folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, 'Top Level Folder'))

        # Find specified type of files
        if self.rb_trdi.isChecked():
            self.get_files(extension='.mmt')

        elif self.rb_sontek.isChecked():
            self.get_files(extension='r.mat')
            self.clean_sontek_files()

        elif self.rb_qrev.isChecked():
            self.get_files(extension='qrev.mat')

        # Enable save button
        self.pb_save.setEnabled(True)

    def get_files(self, extension):
        """Finds all files in the folder and subfolders with the specified extension.

        extension: str
            File extension of files to find.
        """

        for folder_path, folder_names, files in os.walk(self.top_folder):
            for filename in files:
                if filename.lower().endswith(extension):
                    self.files.append(os.path.join(folder_path, filename))

    def clean_sontek_files(self):
        """Removes moving-bed test files from the list."""

        file_list = copy.deepcopy(self.files)
        for filename in file_list:
            if os.path.basename(filename).startswith('Smba'):
                self.files.remove(filename)
            elif os.path.basename(filename).startswith('Loop'):
                self.files.remove(filename)

    def save_list(self):
        """Saves full path and file name to ASCII file, on file per line."""

        # Get output file name from user
        out_file = QtWidgets.QFileDialog.getSaveFileName(self, caption='Save File',
                                                         directory=os.path.join(self.top_folder, 'files.txt'),
                                                         filter='*.txt')[0]
        # Write file
        with open(out_file, 'w') as out_f:
            for file in self.files:
                out_f.write('%s\n' % file)

    def process(self):
        """Process files using QRevPy"""

        # Get name of file to process from user
        fullname = QtWidgets.QFileDialog.getOpenFileNames(
            self, self.tr('Open File'))[0][0]

        # Read entire file of path names, one per line
        with open(fullname, 'r') as f:
            lines = [line.rstrip() for line in f]

        # Process TRDI files
        if lines[0].endswith('.mmt'):
            for line in lines[self.sb_skip.value()::]:
                path, name = os.path.split(line)
                self.txt_status.setText('Processing ' + name)
                QtCore.QCoreApplication.processEvents()
                meas = Measurement(in_file=line, source='TRDI')
                self.append_result(path=path, meas=meas)

        # Process QRev files
        elif lines[0].endswith('QRev.mat'):
            for line in lines[self.sb_skip.value()::]:
                path, name = os.path.split(line)
                self.txt_status.setText('Processing ' + name)
                QtCore.QCoreApplication.processEvents()
                meas = Measurement(in_file=line, source='QRev')
                self.append_result(path=path, meas=meas)

        # Process SonTek files
        elif lines[0].endswith('r.mat'):
            # Find unique paths to identify files that belong to a measurement
            paths = []
            for line in lines[self.sb_skip.value()::]:
                path, name = os.path.split(line)
                paths.append(path)
            unique_path = set(paths)

            # Process each measurement
            for path in unique_path:
                measurement_files = [line for line in lines if path in line]
                folder_name = os.path.basename(path)
                self.txt_status.setText('Processing ' + folder_name)
                QtCore.QCoreApplication.processEvents()
                meas = Measurement(in_file=measurement_files, source='SonTek')
                self.append_result(path=folder_name, meas=meas)
        self.txt_status.setText('Processing Complete')
        self.save_excel(excel_name=self.QExcelFile.text())

    def append_result(self, path, meas):
        for n in range(len(meas.discharge)):
            summary_list = [path,
                            meas.transects[n].file_name,
                            meas.transects[n].boat_vel.selected,
                            meas.discharge[n].total_uncorrected,
                            meas.discharge[n].total,
                            meas.discharge[n].top,
                            meas.discharge[n].middle,
                            meas.discharge[n].bottom,
                            meas.discharge[n].left,
                            meas.discharge[n].right,
                            meas.extrap_fit.sel_fit[-1].top_method,
                            meas.extrap_fit.sel_fit[-1].bot_method,
                            meas.extrap_fit.sel_fit[-1].exponent]
            self.summary.append(summary_list)

    def save_excel(self, excel_name):

        workbook = xlsxwriter.Workbook(excel_name+'.xlsx')
        worksheet = workbook.add_worksheet()
        column_labels = ['Path', 'File', 'Ref', 'Q_unc', 'Q_corr', 'Q_top', 'Q_middle', 'Q_bottom', 'Q_left', 'Q_right',
                         'Top_method', 'Bot_method', 'Exp']
        # Write column labels
        for col, label in enumerate(column_labels):
            worksheet.write(0, col, label)

        # Write data
        for row, sublist in enumerate(self.summary):
            for col, value in enumerate(sublist):
                try:
                    worksheet.write(row + 1, col, value)
                except:
                    worksheet.write(row + 1, col, '')

        # Close workbook
        workbook.close()

app = QtWidgets.QApplication(sys.argv)
window = BatchDialog()
window.show()
app.exec_()
