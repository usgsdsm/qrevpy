

**QRevPy** is a Python port of the Matlab code QRev developed by the USGS to to compute the discharge from a moving-boat ADCP measurement using data collected with any of the Teledyne RD Instrument (TRDI) or SonTek bottom tracking ADCPs. QRev improves the consistency and efficiency of processing streamflow measurements by providing:

* Automated data quality checks with feedback to the user
* Automated data filtering
* Automated application of extrap, LC, and SMBA algorithms
* Consistent processing algorithms independent of the ADCP used to collect the data
* Improved handing of invalid data
* An estimated uncertainty to help guide the user in rating the measurement


**For a full description and instructions on the use of QRev for Matlab click** **[HERE](https://hydroacoustics.usgs.gov/movingboat/QRev.shtml)** **to view the QRev page on the USGS Office of Surface Water Hydroacoustics website.**

***

# Development
**QRevPy** is an active project and is undergoing continuous development. Versions available on this repository may be versions that are being tested. Currently the computational engine is functional but the user interface consist of only a couple of buttons to load and save processed data. Work is ongoing to provide loading of QRev Matlab output and saving an xml file consistent with QRev Matlab. Conceptual mockups of the new user interface are shown in the QRev_Python_User_Interface_Design.pptx file in this repository. If you would like to contribute, please use the pull request process to provide new or improved code. 

## Requirements and Dependencies
### Source Code

QRevPy is currently being developed using Python 3.6.6 and makes use the following packages:

numpy==1.15.2
PyQt5==5.10.1
pytest==3.8.2
scipy==1.1.0
statsmodels==0.9.0
utm==0.4.2
xmltodict==0.11.0


## Bugs
Please report all bugs with appropriate instructions and files to reproduce the issue and add this to the issues tracking feature in this repository.

## Disclaimer
This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.

## License

Copyright / License - CC0 1.0: The person who associated a work with this deed has dedicated the work to the public domain by waiving all of his or her rights to the work worldwide under copyright law, including all related and neighboring rights, to the extent allowed by law. You can copy, modify, distribute and perform the work, even for commercial purposes, all without asking permission. 

In no way are the patent or trademark rights of any person affected by CC0, nor are the rights that other persons may have in the work or in how the work is used, such as publicity or privacy rights.

Unless expressly stated otherwise, the person who associated a work with this deed makes no warranties about the work, and disclaims liability for all uses of the work, to the fullest extent permitted by applicable law.

When using or citing the work, you should not imply endorsement by the author or the affirmer.

Publicity or privacy: The use of a work free of known copyright restrictions may be otherwise regulated or limited. The work or its use may be subject to personal data protection laws, publicity, image, or privacy rights that allow a person to control how their voice, image or likeness is used, or other restrictions or limitations under applicable law.

Endorsement: In some jurisdictions, wrongfully implying that an author, publisher or anyone else endorses your use of a work may be unlawful.

3rd party code is covered by the copyright and license associated with those codes.

## Author
David S Mueller  
U.S. Geological Survey  
9818 Bluegrass Parkway  
Louisville, KY  
<dmueller@usgs.gov>