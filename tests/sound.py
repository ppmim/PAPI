# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 08:57:20 2013

@author: panic
"""

import sys

from PyQt4.QtGui import QApplication, QMainWindow, QDirModel, QColumnView
from PyQt4.QtGui import QFrame
from PyQt4.QtCore import SIGNAL
from PyQt4.phonon import Phonon

class MainWindow(QMainWindow):

    m_model = QDirModel()
    def __init__(self):
        	QMainWindow.__init__(self)
        	self.m_fileView = QColumnView(self)
        	self.m_media = None
    
        	self.setCentralWidget(self.m_fileView)
        	self.m_fileView.setModel(self.m_model)
        	self.m_fileView.setFrameStyle(QFrame.NoFrame)
    
        	self.connect(self.m_fileView,
        		SIGNAL("updatePreviewWidget(const QModelIndex &)"), self.play)

    def play(self, index):
        	self.delayedInit()
        	self.m_media.setCurrentSource(
        		Phonon.MediaSource(self.m_model.filePath(index)))
        	self.m_media.play()

    def delayedInit(self):
        	if not self.m_media:
        		self.m_media = Phonon.MediaObject(self)
        		audioOutput = Phonon.AudioOutput(Phonon.MusicCategory, self)
        		Phonon.createPath(self.m_media, audioOutput)

def main():
        app = QApplication(sys.argv)
        QApplication.setApplicationName("Phonon Tutorial 2 (Python)")
        mw = MainWindow()
        mw.show()
        sys.exit(app.exec_())

if __name__ == '__main__':
    main()
