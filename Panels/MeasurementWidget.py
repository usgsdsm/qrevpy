from PyQt5.QtWidgets import (QWidget, QGridLayout, QSizePolicy, QTableWidget, 
                             QTableWidgetItem, QAbstractItemView, QLabel,
                             QCheckBox, QPushButton, QVBoxLayout, QScrollArea,
                             QComboBox)
from Panels.MatplotlibWidget import MatplotlibWidget
from PyQt5.Qt import QHBoxLayout
from PyQt5.QtCore import QCoreApplication, Qt, QPropertyAnimation,  QRect, QByteArray, QSize
from functools import partial


class MeasurementWidget(QWidget):
    """description of class"""

    def __init__(self, **kargs):
        super(MeasurementWidget, self).__init__()
       
        self.layout = QGridLayout()
        self.expand = []
        self.labels = []
        self.datatables = []
        self.initUI()
        

    def initUI(self):

        
        self.scroll = QScrollArea()
        #self.scroll.setMaximumHeight(550)
        #self.scroll.setMaximumWidth(750)
        
        

        self.a = QWidget()
        self.aLayout = self.grid = QGridLayout()
        
     
        self.a.setLayout(self.aLayout)
        self.a.setMinimumHeight(1100)
        

   
        #labels
        lbl1 = QLabel('<b style="font-size: 12pt">Measurement Details</b>')
        lbl1.move(0,0)
        lbl1.resize(250, 30)
        self.layout.addWidget(lbl1, 0,0,1,1)
        #self.aLayout.addWidget(lbl1)

     
        
        #---------------------------------------------------------------------------------------------------------
        
        h_headers = ['Parameters','Measurement','000_14-04-14','000_14-04-14','002_14-04-14','003_14-04-14']
        v_headers = ['','','','','','','']
        self.chbx1, self.chbx2, self.chbx3, self.chbx4 = QCheckBox(), QCheckBox(), QCheckBox(),QCheckBox()
        data = [['Use','', self.chbx1,self.chbx2, self.chbx3, self.chbx4],
                [QLabel('<b>Total Q (m3/s)</b>'),QLabel('<b>2</b>'),QLabel('<b>3</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>')],
                ['Top Q (m3/s)','5','6','4','5','66'],
                ['Middle Q (m3/s)','2','3','3','5','6'],
                ['Bottom Q (m3/s)','4','4','3','2','1'],
                ['Left Q (m3/s)','2','3','44','2','22'],
                ['Right Q (m3/s)','2','3','3','333','222']]

        self.datatables.append(self.create_table(310, data, h_headers, v_headers))

        func = partial(self.table_toggle, table_idx = 0)
        label, collapse_widget = self.createTableCollapse("Discharge", func)
        self.labels.append(label)
        self.expand.append(True)
        self.aLayout.addWidget(collapse_widget)
        self.aLayout.addWidget(self.datatables[0])
         #---------------------------------------------------------------------------------------------------------

        self.chbx1, self.chbx2, self.chbx3, self.chbx4 = QCheckBox(), QCheckBox(), QCheckBox(),QCheckBox()
        data = [['Use','', self.chbx1,self.chbx2, self.chbx3, self.chbx4],
                [QLabel('<b>Total Q (m3/s)</b>'),QLabel('<b>2</b>'),QLabel('<b>3</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>')],
                ['Top Q (m3/s)','5','6','4','5','66'],
                ['Middle Q (m3/s)','2','3','3','5','6'],
                ['Bottom Q (m3/s)','4','4','3','2','1'],
                ['Left Q (m3/s)','2','3','44','2','22'],
                ['Right Q (m3/s)','2','3','3','333','222']]

        self.datatables.append(self.create_table(310, data, h_headers, v_headers))

        func = partial(self.table_toggle, table_idx=1)
        label,collapse_widget = self.createTableCollapse("The Sequel", func)
        self.labels.append(label)
        self.expand.append(True)
        self.aLayout.addWidget(collapse_widget)
        self.aLayout.addWidget(self.datatables[1])
         #---------------------------------------------------------------------------------------------------------

        self.chbx1, self.chbx2, self.chbx3, self.chbx4 = QCheckBox(), QCheckBox(), QCheckBox(),QCheckBox()
        data = [['Use','', self.chbx1,self.chbx2, self.chbx3, self.chbx4],
                [QLabel('<b>Total Q (m3/s)</b>'),QLabel('<b>2</b>'),QLabel('<b>3</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>')],
                ['Top Q (m3/s)','5','6','4','5','66'],
                ['Middle Q (m3/s)','2','3','3','5','6'],
                ['Bottom Q (m3/s)','4','4','3','2','1'],
                ['Left Q (m3/s)','2','3','44','2','22'],
                ['Right Q (m3/s)','2','3','3','333','222']]

        self.datatables.append(self.create_table(310, data, h_headers, v_headers))

        func = partial(self.table_toggle, table_idx=2)
        label,collapse_widget = self.createTableCollapse("The Triquel", func)
        self.labels.append(label)
        self.expand.append(True)
        self.aLayout.addWidget(collapse_widget)
        self.aLayout.addWidget(self.datatables[2])

        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

        self.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
       

        #self.layout.addWidget(self.scroll)
        #layout = new QVBoxLayout(central);
        #q = QSize()
        #q.setHeight(1000)
        #q.setWidth(700)
        #a.setMaximumSize(q)
        self.scroll.setWidget(self.a)
        self.scroll.setWidgetResizable(True)
       
        
        self.layout.addWidget(self.scroll,1,0,1,1)

        self.scroll2 = QScrollArea()
        self.scroll2.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scroll2.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.chbx1, self.chbx2, self.chbx3, self.chbx4 = QCheckBox(), QCheckBox(), QCheckBox(),QCheckBox()
        data = [['Use','', self.chbx1,self.chbx2, self.chbx3, self.chbx4],
                [QLabel('<b>Total Q (m3/s)</b>'),QLabel('<b>2</b>'),QLabel('<b>3</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>'),QLabel('<b>9</b>')],
                ['Top Q (m3/s)','5','6','4','5','66'],
                ['Middle Q (m3/s)','2','3','3','5','6'],
                ['Bottom Q (m3/s)','4','4','3','2','1'],
                ['Left Q (m3/s)','2','3','44','2','22'],
                ['Right Q (m3/s)','2','3','3','333','222']]

        
        self.scroll2.setWidgetResizable(True)
        #labels

        lbl2 = QLabel('<b style="font-size: 12pt">Messages</b>')
        lbl2.resize(250, 30)
        self.layout.addWidget(lbl2,2,0,1,1)
       
        #self.aLayout.addWidget(lbl1)
        
        self.b = QWidget()
        self.bLayout = QVBoxLayout()
        self.b.setLayout(self.bLayout)
        self.bLayout.addWidget(self.create_table(200, data, h_headers, v_headers))
        self.scroll2.setWidget(self.b)
        self.layout.addWidget(self.scroll2,3,0,1,1)


        self.rightPane = QWidget()
        self.rightLayout = QVBoxLayout()
        self.rightPane.setLayout(self.rightLayout)

        #labels
        lbl3 = QLabel('<b style="font-size: 12pt">Measurement Quality Assessment</b>')
        lbl3.resize(250, 30)
        self.rightLayout.addWidget(lbl3)


        h_headers = ['','COV %','','% Q' ]
        v_headers = ['','','']
        data = [
                ['Q:','11.92','Left/Right Edge:','21.81 / 4.03'],
                ['Width:','8.73','Invalid Cells:','0.64'],
                ['Area:','11.32','Invalid Ens:','0.08',]]

        self.datatables.append(self.create_table(110, data, h_headers, v_headers))
        self.datatables[3].setColumnWidth(0, 50)
        self.datatables[3].setColumnWidth(1, 60)
        self.datatables[3].setColumnWidth(2,100)

        for x in range(len(data)):
            self.datatables[3].setRowHeight(x,25)
      
        self.rightLayout.addWidget(self.datatables[3])

        h_headers = ['Parameter','Automatic','User']
        v_headers = ['','','','','','','']
        data = [['Random 95% Uncertainty', '19.0', ''],
                 ['Invalid Data 95% Uncertainty', '0.1', ''],	
                ['Edge Q 95% Uncertainty',	'7.8', ''],
                ['Extrapolation 95% Uncertainty', '0.6', '' ],
                ['Moving-Bed 95% Uncertainty', '1.5', ''],
                ['Systematic 68% Uncertainty', '1.5', ''],
                ['Estimated 95% Uncertainty', '20.8','20.8']]
        

        self.datatables.append(self.create_table(210, data, h_headers, v_headers))
        for x in range(len(data)):
            self.datatables[4].setRowHeight(x,25)
        
      
        self.rightLayout.addWidget(self.datatables[4])

      

        lbl4 = QLabel('<b style="font-size: 12pt">User Rating</b>')
        lbl4.resize(250, 30)
        self.rightLayout.addWidget(lbl4)

        self.combo = QComboBox()
        self.combo.addItem('Not Rated')
        self.combo.addItem('Excellent (<2%)')
        self.combo.addItem('Good (2-5%)')
        self.combo.addItem('Fair (5-8%)')
        self.combo.addItem('Poor (>8%)')
        self.rightLayout.addWidget(self.combo)
         
        

        lbl5 = QLabel('<b style="font-size: 12pt">Profile Extrapolation</b>')
        lbl5.resize(250, 30)
        self.rightLayout.addWidget(lbl5)
       
        x=[0,10,100]
        y=[3,4,5]

        self.mplwidget = MatplotlibWidget()

        
       
        self.mplwidget.setObjectName("mplwidget")
        
        self.mplwidget.plotDataPoints(x,y)
       
        
        self.rightLayout.addWidget(self.mplwidget);

        self.rightLayout.setStretch(0,1)
        self.rightLayout.setStretch(1,3)
        self.rightLayout.setStretch(2,5)
        self.rightLayout.setStretch(3,1)
        self.rightLayout.setStretch(4,1)
        self.rightLayout.setStretch(5,1)
        self.rightLayout.setStretch(6,8)

        self.layout.addWidget(self.rightPane,0,1,4,1)
       
        
        #q = QSize()
        #q.setHeight(1000)
        #q.setWidth(700)
        #a.setMaximumSize(q)

        #self.scroll.setLayout(self.layout)
        self.layout.setRowStretch(0,1)
        self.layout.setRowStretch(1, 10)
        self.layout.setRowStretch(2, 1)
        self.layout.setColumnStretch(0, 2)
        self.layout.setColumnStretch(1, 1)

        self.setLayout(self.layout)
      
    def createTableCollapse(self, text, func):

        collapseWidget = QWidget()
        hlayout = QHBoxLayout()
        collapseWidget.setLayout(hlayout)
        button = QPushButton('-')
       
       
        hlayout.setSpacing(10)

        button.clicked.connect(func)
        button.setFixedWidth(40)
        
        label = QLabel('<b style="font-size: 10pt">%s</b>' % text)
       
        hlayout.addWidget(button, 1)
        hlayout.addWidget(label, 20)

        return (button, collapseWidget)

    def create_table(self, height, data, h_headers, v_headers = None):
        datatable = QTableWidget()
        datatable.setMaximumHeight(height)
        datatable.setRowCount(len(data))
        datatable.setColumnCount(len(data[0]))
        
        for x in range(len(data)):
            for y in range(len(data[x])):
                if type(data[x][y]) is str:
                    datatable.setItem(x, y, QTableWidgetItem(data[x][y]))
                else:
                    datatable.setCellWidget(x,y,data[x][y])

        datatable.setHorizontalHeaderLabels(h_headers)
        if v_headers is not None:
            datatable.setVerticalHeaderLabels(v_headers)

        datatable.setAlternatingRowColors(True)
        datatable.setSortingEnabled(False)
        datatable.setEditTriggers(QAbstractItemView.NoEditTriggers)
        datatable.setStyleSheet("QHeaderView::section { background-color: #e5e5e5; font-size: 8pt; font-weight: bold} "
        "QTableWidget { font-size: 8pt; } QHeaderView::section:horizontal { height: 30px}")
        datatable.horizontalHeader().setStretchLastSection(True);

        return datatable




    def table_toggle(self, table_idx):

        

        dt = self.datatables[table_idx]

        def resize_layout():
            dt.setMaximumHeight(0)

        if self.expand[table_idx] == True:
            animation = QPropertyAnimation(self)
            animation.setPropertyName(QByteArray().append('size'))
            animation.setTargetObject(dt)
            animation.setDuration(250)

            ogWidth = dt.width()
            animation.setStartValue(QSize(ogWidth,310))
            animation.setEndValue(QSize(ogWidth,0))
            animation.finished.connect(resize_layout)
            animation.start()
           
            self.labels[table_idx].setText('+')
            
            self.a.updateGeometry()
            self.a.setMinimumHeight(self.a.minimumHeight() - 320)
            
           
        else:
            animation = QPropertyAnimation(self)
            animation.setPropertyName(QByteArray().append('size'))
            animation.setTargetObject(dt)
            animation.setDuration(250)

            dt.setMaximumHeight(310)
            ogWidth = dt.width()
            animation.setStartValue(QSize(ogWidth,0))
            animation.setEndValue(QSize(ogWidth,310))
            animation.start()

            self.labels[table_idx].setText('-')
           
            self.a.updateGeometry()
            self.a.setMinimumHeight(self.a.minimumHeight() + 320)



        self.expand[table_idx] = not self.expand[table_idx]



