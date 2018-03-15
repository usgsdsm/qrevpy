import scipy.io as sio


class MatSonTek(object):
    """Read SonTek Matlab files and returns a dictionary of mat_struct.
     Any data in English units are converted to SI units.

     Atributes
     ---------
     mat_data: dictionary
        This is a dictionary of mat_struct, which can be accessed as mat_data[key].struct_name.sub_name.
    convert2correct: method
        Converts all English data to metric units.
    """

    def __init__(self, fullname):
        """Initializes the object, reads the Matlab file, and converts all English units to metric.

        Parameters
        ----------
        fullname: str
            String contain both the path and filename.
        """

        # Read Matlab file
        mat_data = sio.loadmat(fullname, struct_as_record=False, squeeze_me=True)

        # Convert data to SI units if in English units
        if mat_data['BottomTrack'].Units.BT_Depth == 'ft':
            self.convert2metric(mat_data)

        vars(self).update(mat_data)


    def convert2metric(self, mat_data):
        """Converts all data in English units to metric units.
        """

        data2correct = ['BottomTrack', 'GPS', 'Setup', 'Summary', 'System', 'WaterTrack']
        for item in data2correct:
            data = mat_data[item]
            units = data.Units
            names = units._fieldnames
            for name in names:
                if getattr(units, name) == 'ft':
                    setattr(data, name, getattr(data, name) * 0.3048)
                    setattr(units, name, 'm')
                elif getattr(units, name) == 'ft/s':
                    setattr(data, name, getattr(data, name) * 0.3048)
                    setattr(units, name, 'm/s')
                elif getattr(units, name) == 'degF':
                    setattr(data, name, (getattr(data, name)-32) * (5.0/9.0))
                    setattr(units, name, 'degC')
                elif getattr(units, name) == 'cfs':
                    setattr(data, name, getattr(data, name) * (0.3048**3))
                    setattr(units, name, 'm3/s')
                elif getattr(units, name) == 'ft2':
                    setattr(data, name, getattr(data, name) * (0.3048 ** 2))
                    setattr(units, name, 'm2')
