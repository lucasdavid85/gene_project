#!/usr/bin/env python3

#requirements : 
#pip install pyqt6
#pip install Biopython





#Library importation
import os.path
import networkx as nx 
from Bio import pairwise2
import matplotlib.pyplot as plt
import sys
import subprocess
import os.path
from pathlib import Path
from PyQt6 import QtWidgets
from pyqtgraph import PlotWidget,plot
import pyqtgraph as pg


#from PyQt6.QtWidgets import (
#    QApplication,
#    QDialog,
#    QDialogButtonBox,
#    QFormLayout,
#    QLineEdit,
#    QVBoxLayout,
#    QLabel,
#    QPushButton,
#    QWidget,
#)


#class Window(QDialog):
#
#    def __init__(self):
#        super().__init__(parent=None)
#        self.setWindowTitle("QDialog")
#        dialogLayout = QVBoxLayout()
#        formLayout = QFormLayout()
##        formLayout.addRow("Job:", QLineEdit())
##        formLayout.addRow("Hobbies:", QLineEdit())
#        dialogLayout.addLayout(formLayout)
#        buttons = QDialogButtonBox()
#        buttons.setStandardButtons(
#            QDialogButtonBox.StandardButton.Cancel
#            | QDialogButtonBox.StandardButton.Ok
#        )
#        dialogLayout.addWidget(buttons)
#        self.setLayout(dialogLayout)
#
#
#
#if __name__ == "__main__":
#
#    app = QApplication([])
#
#    window = Window()
#
#    window.show()
#
#    sys.exit(app.exec())

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):



        #Get fasta file name
        fasta_name=sys.argv[1]
        file_exist=os.path.exists(fasta_name)
        
        if(len(sys.argv)!=2):
            print("Please use this type of syntax : python3 M2BB_Project_Lucas_David.py <name_of_fasta_file>")
            exit(0)


        if(file_exist):
            print("The file exist\n")
            fasta_dic=(fasta_to_dic(fasta_name))

            # Ask the name of the folder and files  
            filename=input("Please enter the filename for the graph without extentions\n")

            #Make the pairwise
            make_pairwise(fasta_dic,filename)
        else:
            print("The file does not exist or is empty or has been wrong writing, please enter a new filename")
            exit(0)

        #Ask if the user want to open the graph into Gephi software
        open_gephi_yn=input("Do you want to open Gephi to visualize the graph ? y for yes, n for no \n")

        if(open_gephi_yn=="y"):
            gephi_filename="{}/{}.gexf".format(filename,filename)
            home_path=Path.home()
            path_to_gephi="{}/gephi-0.9.7/bin/gephi".format(home_path)
            result=subprocess.run([path_to_gephi, "-o",gephi_filename])






#        super(MainWindow, self).__init__(*args, **kwargs)
#
#        self.graphWidget = pg.PlotWidget()
#        self.setCentralWidget(self.graphWidget)
#
#        hour = [1,2,3,4,5,6,7,8,9,10]
#        temperature = [30,32,34,32,33,31,29,32,35,45]
#
#        # plot data: x, y values
#        self.graphWidget.plot(hour, temperature)


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()









#def greet():
#
#
#
#    if msgLabel.text():
#        msgLabel.setText("")
#    else:
#        msgLabel.setText("Hello, World!")
#
#app = QApplication([])
#window = QWidget()
#window.setWindowTitle("Signals and slots")
#layout = QVBoxLayout()
#button = QPushButton("Greet")
#button.clicked.connect(greet)
#layout.addWidget(button)
#msgLabel = QLabel("")
#layout.addWidget(msgLabel)
#window.setLayout(layout)
#window.show()
#sys.exit(app.exec())











#verifying the number of arguments



def fasta_to_dic(fasta_name): 
    '''Take fasta name and transform it from a list to a dictionnary'''
    fasta_file=open(fasta_name,'r')
    fasta_dic={}
    for line in fasta_file:
        if line[0]=='>':
            header=line[1:-1]
            header=header[:25]

            fasta_dic[header]=''
        else:
            fasta_dic[header]+=line[:-1]
    fasta_file.close()
    return(fasta_dic)

def make_pairwise(fasta_dic,filename):
    '''Take the fasta dictionnary created and make the pairwise beetween all the sequences in the fasta file This definition retun the score of all pairwise alignments and create a graph, saved in png and gexf'''
    
    #Get the list of the keys
    fasta_list_keys=list(fasta_dic.keys())

    #get the list of the values
    fasta_list_values=list(fasta_dic.values())
    
    alignments=[]

    G= nx.Graph()
    for i in range(len(fasta_list_keys)):
        G.add_node(fasta_list_keys[i])

    G.add_nodes_from(fasta_list_keys)
    compt=0
    for i in range(0,len(fasta_list_values)):
        #print(fasta_list_keys[i])
        for j in range(i+1,len(fasta_list_values)):
            alignments.append(pairwise2.align.globalxx(fasta_list_values[i],fasta_list_values[j],score_only=True))
            G.add_edge(fasta_list_keys[i],fasta_list_keys[j],weight=alignments[compt])
            compt+=1

    #Drawing the graph
    pos=nx.circular_layout(G)
    labels=nx.get_edge_attributes(G,'weight')
    nx.draw(G,pos,with_labels=True,node_size=15, font_size = 10,node_color="blue",width=0.5)
    esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.5]
    nx.draw_networkx_edges(G, pos, edgelist=esmall, width=10, alpha=0.5, edge_color="b")
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
    plt.axis("on")

    #Create folder
    subprocess.run(["mkdir",filename])

    #Saving all files
    nx.write_gexf(G,f'{filename}/{filename}.gexf')
    print("Saving the gefx file to {}/{}.gexf...".format(filename,filename))
    plt.savefig(f"{filename}/{filename}.png")
    print("Saving the gefx file to {}/{}.png...".format(filename,filename))
    plt.show()
    return alignments


##### Starting the program ######
#verify if the fasta_file exist, if is not empty and if it's the right extention











