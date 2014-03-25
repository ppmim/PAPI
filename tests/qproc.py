import sys
from qt import qApp, QApplication, QGridLayout, QLineEdit, QProcess, \
               QPushButton, QString, QStringList, QTextBrowser, QTimer, \
               QWidget, SIGNAL, SLOT

class Window(QWidget):


    def __init__(self):
    
        QWidget.__init__(self)
        
        self.textBrowser = QTextBrowser(self)
        self.textBrowser.setTextFormat(QTextBrowser.LogText)
        self.lineEdit = QLineEdit(self)
        self.startButton = QPushButton(self.tr("Start"), self)
        self.stopButton = QPushButton(self.tr("Stop"), self)
        self.stopButton.setEnabled(False)

        self.connect(self.lineEdit, SIGNAL("returnPressed()"), self.startCommand)
        self.connect(self.startButton, SIGNAL("clicked()"), self.startCommand)
        self.connect(self.stopButton, SIGNAL("clicked()"), self.stopCommand)
        layout = QGridLayout(self, 2, 3)
        layout.setSpacing(8)
        layout.addMultiCellWidget(self.textBrowser, 0, 0, 0, 2)
        layout.addWidget(self.lineEdit, 1, 0)
        layout.addWidget(self.startButton, 1, 1)
        layout.addWidget(self.stopButton, 1, 2)
        self.process = QProcess()
        self.connect(self.process, SIGNAL("readyReadStdout()"), self.readOutput)
        self.connect(self.process, SIGNAL("readyReadStderr()"), self.readErrors)
        self.connect(self.process, SIGNAL("processExited()"), self.resetButtons)

    def startCommand(self):
        self.process.setArguments(QStringList.split(" ", self.lineEdit.text()))
        self.process.closeStdin()
        self.startButton.setEnabled(False)
        self.stopButton.setEnabled(True)
        self.textBrowser.clear()

        if not self.process.start():
            self.textBrowser.setText(
                QString("*** Failed to run %1 ***").arg(self.lineEdit.text())
                )
            self.resetButtons()
            return

    def stopCommand(self):
        self.resetButtons()
        self.process.tryTerminate()
        QTimer.singleShot(5000, self.process, SLOT("kill()"))

    def readOutput(self):
    
        self.textBrowser.append(QString(self.process.readStdout()))
    
    def readErrors(self):
    
        self.textBrowser.append("error: " + QString(self.process.readLineStderr()))

    def resetButtons(self):
        self.startButton.setEnabled(True)
        self.stopButton.setEnabled(False)

if __name__ == "__main__":

    app = QApplication(sys.argv)
    window = Window()
    app.setMainWidget(window)
    window.show()

    sys.exit(app.exec_loop())