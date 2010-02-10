
import sys
from qt import QWidget, QProcess, \
               QString, QStringList, QTimer, \
               SIGNAL, SLOT

#Log
import misc.paLog
from misc.paLog import log
import reduce # for TaskInfo
  
class RunQtProcess(QWidget):


    def __init__(self, commandToRun, guiOutWindow, taskInfoList , output=None):
    
        QWidget.__init__(self)
        
        self._command = commandToRun
        self._outWindow = guiOutWindow
        self._task_info_list = taskInfoList
        self._cmd_output = output  # Output (file) that the command execution should generate if success
        
        self.exit = 0             # This variable show if the last execution end sucessful (0) or not (1)
        
        self.taskInfo = reduce.TaskInfo()
        
        try:
            self.process = QProcess()
            self.connect(self.process, SIGNAL("readyReadStdout()"), self.readOutput)
            self.connect(self.process, SIGNAL("readyReadStderr()"), self.readErrors)
            self.connect(self.process, SIGNAL("processExited()"), self.exitFunc)
        except:
            raise
        
    def startCommand(self):
        print "RunQtProcess-Command=", self._command
        
        self.process.setArguments(QStringList.split(" ", self._command))
        self.process.closeStdin()

        if not self.process.start():
            self._outWindow.setText(
                QString("*** Failed to run %1 ***").arg(self._command)
                )
            self.exit = 1
            self.exitFunc()    
            return
        
        print "End startCommand"

    def stopCommand(self):
        self.process.tryTerminate()
        QTimer.singleShot(5000, self.process, SLOT("kill()"))

    def readOutput(self):
    
        """
        TO IMPROVE !!!!
        """
        err = str(self.process.readStdout())
        
        if (err.count('error') or err.count('ERROR') or err.count('Segmentation fault') or err.count("command not found")
          or err.count("No such file or directory") or err.count('No match') or err.count("Failed") or err.count("fail")):
            log.error( "An error happened while running command --> %s \n", err)
            self.exit = 1
        else:
          pass 
          #log.info("readOutput: no error detected!")
            
        self._outWindow.append(QString("STDOUT>>> %1").arg(err))
    
    def readErrors(self):
        """
        TO IMPROVE !!!!
        """
        err = str(self.process.readLineStderr())
        
        if (err.count('error') or err.count('ERROR') or err.count('Segmentation fault') or err.count("command not found")
          or err.count("No such file or directory") or err.count('No match') or err.count("Failed") or err.count("fail") or err.count("cannot")):
            log.error( "An error happened while running command --> %s \n", err)
            self.exit = 1
        else:
          pass 
          #log.info("readErrors: NO error detected!")
        
        self._outWindow.append( QString("STDERROR>>> %1").arg(err) )
    
    def exitFunc(self):
        
        self.taskInfo._name =  "RunQtProcess : " + self._command             
        self.taskInfo._exit_status =  self.exit      
        self.taskInfo._return = self._cmd_output           
        self.taskInfo._exc = None
        
        self._task_info_list.append(self.taskInfo)
        
    def normalEnd(self):
              
        return self._exit
    
