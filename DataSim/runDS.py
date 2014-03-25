#!/usr/bin/env python

################################################################################
#
# runDS (run Data Simulator for PANIC pipeline)
#
# runDS.py
#
# Last update 11/Sep/2008
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################


from qt import *

from datasimulator import *
import sys
import os


class MyDS(DataSimulator):
    def source_directory(self):
        dir=QFileDialog.getExistingDirectory("/disk-a/caha/", self,
                                             "get existing directory", "Choose a directory",True )
        self.lineEdit_sourceD.setText(dir)

    def output_directory(self):
        dir=QFileDialog.getExistingDirectory("/disk-a/caha/", self,
                                             "get existing directory", "Choose a directory",True )
        self.lineEdit_outputD.setText(dir)  


    def type_changed(self):

        if self.comboBox_type.currentText()=="Custom":
            self.lineEdit_type.setEnabled(True)
        else:
            self.lineEdit_type.setEnabled(False)

    def run(self):

        print "Running..."
        #SOURCE
        orig=self.lineEdit_sourceD.text()
        #DEST
        dest=self.lineEdit_outputD.text()
        #TYPE
        if self.comboBox_type.currentText()=="Custom":
            type=self.lineEdit_type.text()
        elif self.comboBox_type.currentText()=="All":
            type='all'
        else:
            type=self.comboBox_type.currentText()

        #DELAY
        delay=self.lineEdit_delay.text().toInt()
        #MEF
        if self.checkBox_MEF.isChecked():
            mef="--mef"
        else:
            mef=""
        #TEST
        if self.checkBox_test.isChecked():
            test="--test"
        else:
            test=""
        
        os.system(("./datsimu.py --source=%s --des=%s --type=%s --delay=%d  %s %s" % (orig,dest,type,delay[0],mef,test)))
        

    def stop(self):
        print "Stopping ..."

    def exit(self):
        sys.exit(0)

        
 ################################################################################
if __name__ == "__main__":
    app = QApplication(sys.argv)
    f = MyDS()
    f.show()
    app.setMainWidget(f)
    app.exec_loop()
