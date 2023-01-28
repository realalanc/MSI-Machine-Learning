from PyQt5.QtWidgets import QDialog
from PyQt5 import QtWidgets
from PyQt5 import uic
import sys

class Mainui(QDialog):
    def __init__(self):
        super(Mainui,self).__init__()
        uic.loadUi("causalfind.ui",self)



def initshow():
    app=QtWidgets.QApplication(sys.argv)
    mui=Mainui()
    mui.show()
    sys.exit(app.exec())

if __name__=="__main__":
    initshow()
