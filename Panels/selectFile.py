import os
from PyQt5 import QtWidgets
from Classes.stickysettings import StickySettings as SSet
import wSelectFile
import datetime as datetime

class OpenMeasurementDialog(QtWidgets.QDialog, wSelectFile.Ui_selectFile):
    """Dialog to allow users to select measurement files for processing.

    Parameters
    ----------
    wSelectFile.Ui_selectFile : QDialog
        Dialog window with options for users

    Attributes
    ----------
    settings: dict
        Dictionary used to store user defined settings.
    fullName: list
        Full name of files including path.
    fileName: list
        List of one or more fileNames to be processed.
    pathName: str
        Path to folder containing files.
    type: str
        Type of file (SonTek, TRDI, QRev).
    checked: bool
        Switch for TRDI files (True: load only checked, False: load all).
    pbSonTek: QWidgets.QPushButton
        Allows selection of SonTek Matlab transect files.
    pbTRDI: QWidgets.QPushButton
        Allows selection of TRDI mmt file.
    cbTRDI: QWidgets.QCheckBox
        When checked only the checked transects in the mmt file are loaded.
    pbQRev: QWidgets.QPushButton
        Allows selection of QRev file.
    pbCancel: QWidgets.QPushButton
        Closes dialog with no file selection.
    pbHelp: QWidgets.QPushButton
        Opens help file.
    """

    def __init__(self, parent=None):
        """Initializes settings and connections.

        Parameters
        ----------
        parent
            Identifies parent GUI.
        """

        super(OpenMeasurementDialog, self).__init__(parent)
        self.setupUi(self)

        # Create settings object which contains the default folder
        self.settings = SSet(parent.settingsFile)

        # Create connections for buttons
        self.pbSonTek.clicked.connect(self.selectSonTek)
        self.pbTRDI.clicked.connect(self.selectTRDI)
        self.pbQRev.clicked.connect(self.selectQRev)
        self.pbCancel.clicked.connect(self.cancel)

        # Initialize parameters
        self.fullName = []
        self.fileName = []
        self.pathName = []
        self.type = ''
        self.checked = False

    def defaultFolder(self):
        """Returns default folder.

        Returns the folder stored in settings or if no folder is stored, then the current
        working folder is returned.
        """
        try:
            folder = self.settings.get('Folder')
            if not folder:
                folder = os.getcwd()
        except KeyError:
            self.settings.new('Folder', os.getcwd())
            folder = self.settings.get('Folder')
        return folder

    def processNames(self):
        """Parses fullnames into filenames and pathnames and sets default folder.
        """
        # Parse filenames and pathname from fullName
        if isinstance(self.fullName, str):
            self.pathName, self.fileName = os.path.split(self.fullName)
        else:
            for file in self.fullName:
                self.pathName, fileTemp = os.path.split(file)
                self.fileName.append(fileTemp)

        # Update the folder setting
        self.settings.set('Folder', self.pathName)

    def selectSonTek(self):
        """Get filenames and pathname for SonTek Matlab transect files

        Allows the user to select one or more SonTek Matlab transect files for
        processing. The selected folder becomes the default folder for subsequent
        selectFile requests.
        """

        # Get the current folder setting.
        folder = self.defaultFolder()

        # Get the full names (path + file) of the selected files
        self.fullName = QtWidgets.QFileDialog.getOpenFileNames(
                    self, self.tr('Open File'), folder,
                    self.tr('SonTek Matlab File (*.mat)'))[0]

        # Initialize parameters
        self.type = 'SonTek'
        self.checked = False

        # Process fullName if selection was made
        if self.fullName:
            self.processNames()
        self.close()


        # Read Matlab file
        # mat_contents = sio.loadmat(fileName, struct_as_record=False, squeeze_me=True)

    def selectTRDI(self):
        """Get filenames and pathname for TRDI mmt file

        Allows the user to select a TRDI mmt file for processing.
        The selected folder becomes the default folder for subsequent
        selectFile requests.
        """

        # Get the current folder setting.
        folder = self.defaultFolder()

        # Get the full names (path + file) of the selected files
        self.fullName = QtWidgets.QFileDialog.getOpenFileNames(
                    self, self.tr('Open File'), folder,
                    self.tr('TRDI mmt File (*.mmt)'))[0]

        # Initialize parameters
        self.type = 'TRDI'
        self.checked = self.cbTRDI.isChecked()

        # Process fullName if selection was made
        if self.fullName:
            self.processNames()
        self.close()

    def selectQRev(self):
        """Get filename and pathname of QRev file.

                Allows the user to select a QRev file for viewing or reprocessing.
                The selected folder becomes the default folder for subsequent
                selectFile requests.
                """

        # Get the current folder setting.
        folder = self.defaultFolder()

        # Get the full names (path + file) of the selected file
        self.fullName = QtWidgets.QFileDialog.getOpenFileName(
            self, self.tr('Open File'), folder,
            self.tr('QRev File (*_QRev.mat)'))[0]

        # Initialize parameters
        self.type = 'QRev'
        self.checked = False

        # Process fullName if selection was made
        if self.fullName:
            self.processNames()
        self.close()

    def cancel(self):
        """Close dialog.
        """
        self.close()


class SaveMeasurementDialog(QtWidgets.QDialog):
    """Dialog to allow users to select measurement files for processing.

        Parameters
        ----------
        wSelectFile.Ui_selectFile : QDialog
            Dialog window with options for users

        Attributes
        ----------
        full_Name: str
            Filename with path to save file.
    """

    def __init__(self, parent=None):
        """Initializes settings and connections.

        Parameters
        ----------
        parent
            Identifies parent GUI.
        """
        super(SaveMeasurementDialog, self).__init__(parent)
        # self.setupUi(self)

        # Create settings object which contains the default folder
        settings = SSet('TestSettings')

        # Get the current folder setting.
        folder = self.defaultFolder(settings)

        # Create default file name
        file_name = os.path.join(folder, datetime.datetime.today().strftime('%Y%m%d_%H%M%S_py_QRev.mat'))
        # Get the full names (path + file) of the selected file
        self.full_Name = QtWidgets.QFileDialog.getSaveFileName(
            self, self.tr('Save File'), file_name,
            self.tr('QRev File (*_QRev.mat)'))[0]

    @staticmethod
    def defaultFolder(settings):
        """Returns default folder.

        Returns the folder stored in settings or if no folder is stored, then the current
        working folder is returned.
        """
        try:
            folder = settings.get('Folder')
            if not folder:
                folder = os.getcwd()
        except KeyError:
            settings.new('Folder', os.getcwd())
            folder = settings.get('Folder')
        return folder