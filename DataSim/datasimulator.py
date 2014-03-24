# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'datasimulator.ui'
#
# Created: Thu Sep 11 14:11:13 2008
#      by: The PyQt User Interface Compiler (pyuic) 3.14.1
#
# WARNING! All changes made in this file will be lost!


from qt import *


class DataSimulator(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("DataSimulator")



        self.textLabel2 = QLabel(self,"textLabel2")
        self.textLabel2.setGeometry(QRect(100,112,40,20))

        self.textLabel3 = QLabel(self,"textLabel3")
        self.textLabel3.setGeometry(QRect(80,150,64,20))

        self.pushButton_Stop = QPushButton(self,"pushButton_Stop")
        self.pushButton_Stop.setGeometry(QRect(260,280,108,23))
        self.pushButton_Stop.setPaletteBackgroundColor(QColor(244,65,25))

        self.pushButton_outputD = QPushButton(self,"pushButton_outputD")
        self.pushButton_outputD.setGeometry(QRect(20,70,108,23))

        self.pushButton_sourceD = QPushButton(self,"pushButton_sourceD")
        self.pushButton_sourceD.setGeometry(QRect(20,30,108,23))

        self.pushButton_Run = QPushButton(self,"pushButton_Run")
        self.pushButton_Run.setGeometry(QRect(110,280,108,23))
        self.pushButton_Run.setPaletteBackgroundColor(QColor(69,244,34))

        self.pushButton_Exit = QPushButton(self,"pushButton_Exit")
        self.pushButton_Exit.setGeometry(QRect(410,280,108,23))
        self.pushButton_Exit.setPaletteBackgroundColor(QColor(244,236,92))

        self.lineEdit_delay = QLineEdit(self,"lineEdit_delay")
        self.lineEdit_delay.setGeometry(QRect(150,150,50,21))
        self.lineEdit_delay.setPaletteBackgroundColor(QColor(255,255,155))

        self.checkBox_test = QCheckBox(self,"checkBox_test")
        self.checkBox_test.setGeometry(QRect(340,150,84,20))

        self.lineEdit_sourceD = QLineEdit(self,"lineEdit_sourceD")
        self.lineEdit_sourceD.setGeometry(QRect(150,30,320,21))
        self.lineEdit_sourceD.setPaletteBackgroundColor(QColor(255,255,155))

        self.checkBox_MEF = QCheckBox(self,"checkBox_MEF")
        self.checkBox_MEF.setGeometry(QRect(250,150,84,20))

        self.comboBox_type = QComboBox(0,self,"comboBox_type")
        self.comboBox_type.setGeometry(QRect(152,110,78,21))

        self.lineEdit_outputD = QLineEdit(self,"lineEdit_outputD")
        self.lineEdit_outputD.setGeometry(QRect(150,70,320,21))
        self.lineEdit_outputD.setPaletteBackgroundColor(QColor(255,255,155))

        self.lineEdit_type = QLineEdit(self,"lineEdit_type")
        self.lineEdit_type.setEnabled(0)
        self.lineEdit_type.setGeometry(QRect(250,110,89,21))
        self.lineEdit_type.setPaletteBackgroundColor(QColor(255,255,155))

        self.languageChange()

        self.resize(QSize(546,349).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.pushButton_sourceD,SIGNAL("clicked()"),self.source_directory)
        self.connect(self.pushButton_outputD,SIGNAL("clicked()"),self.output_directory)
        self.connect(self.pushButton_Run,SIGNAL("clicked()"),self.run)
        self.connect(self.pushButton_Exit,SIGNAL("clicked()"),self.exit)
        self.connect(self.comboBox_type,SIGNAL("activated(int)"),self.type_changed)
        self.connect(self.pushButton_Stop,SIGNAL("clicked()"),self.stop)


    def languageChange(self):
        self.setCaption(self.__tr("PANIC Data Simulator"))
        self.textLabel2.setText(self.__tr("Type"))
        self.textLabel3.setText(self.__tr("Delay (sec)"))
        self.pushButton_Stop.setText(self.__tr("Stop"))
        self.pushButton_outputD.setText(self.__tr("Output directory ..."))
        self.pushButton_sourceD.setText(self.__tr("Source directory ..."))
        self.pushButton_Run.setText(self.__tr("Run"))
        self.pushButton_Exit.setText(self.__tr("Exit"))
        self.lineEdit_delay.setText(self.__tr("2"))
        self.checkBox_test.setText(self.__tr("Test"))
        self.checkBox_MEF.setText(self.__tr("MEF"))
        self.comboBox_type.clear()
        self.comboBox_type.insertItem(self.__tr("All"))
        self.comboBox_type.insertItem(self.__tr("Dark"))
        self.comboBox_type.insertItem(self.__tr("Flats"))
        self.comboBox_type.insertItem(self.__tr("Sky Flats"))
        self.comboBox_type.insertItem(self.__tr("Dome Flats"))
        self.comboBox_type.insertItem(self.__tr("Science"))
        self.comboBox_type.insertItem(self.__tr("Custom"))


    def source_directory(self):
        print "DataSimulator.source_directory(): Not implemented yet"

    def output_directory(self):
        print "DataSimulator.output_directory(): Not implemented yet"

    def run(self):
        print "DataSimulator.run(): Not implemented yet"

    def stop(self):
        print "DataSimulator.stop(): Not implemented yet"

    def exit(self):
        exit(1)
        

    def type_changed(self):
        print "DataSimulator.type_changed(): Not implemented yet"

    def __tr(self,s,c = None):
        return qApp.translate("DataSimulator",s,c)
