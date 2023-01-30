from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets
from PyQt5 import uic
import sys

class Mainui(QMainWindow):
    def __init__(self):
        super(Mainui,self).__init__()
        uic.loadUi("causalfind.ui",self)
    
    def pb2_click(self):
        og=self.mainfind()

    def load_csv(self,filename):
        import pandas as pd
        self.data=pd.read_csv(filename)

    def load_h5(self,filename):
        import anndata as ad
        data=ad.read_h5ad(filename)

    def mainfind(self):
        import cdt
        import networkx as nx
        import matplotlib.pyplot as plt
        glasso = cdt.independence.graph.Glasso()
        skeleton = glasso.predict(self.datat)
        new_skeleton = cdt.utils.graph.remove_indirect_links(skeleton, alg='aracne')
        model = cdt.causality.graph.PC()
        output_graph=model.predict(self.datat,new_skeleton)
        return output_graph


if __name__=="__main__":
    app=QtWidgets.QApplication(sys.argv)
    mui=Mainui()
    mui.show()
    mui.pushButton_2.clicked.connect(mui.pb2_click)
    app.exec_()
