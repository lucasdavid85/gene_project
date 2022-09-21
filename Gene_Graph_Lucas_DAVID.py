#!/usr/bin/env python3
#Library importation
import sys
#from turtle import filling
from PyQt5.QtWidgets import QApplication, QWidget,QPushButton,QHBoxLayout,QFileDialog
import os.path
import networkx as nx 
from Bio import pairwise2
import matplotlib.pyplot as plt
import sys
import subprocess
import os.path
from pathlib import Path
import signal



#Notice that you have to download gephi in the home path


class Fenetre(QWidget):

    #To visualize the window
    def __init__(self):
        QWidget.__init__(self)
        dialog = QFileDialog()
        self.fname = dialog.getOpenFileName(None, "Import FASTA", "", "FASTA data files (*.fasta)")
        filename=self.fname[0]
        file=filename[:-6]
        head, tail =os.path.split(file)
        fasta_file=open(filename,'r')
        fasta_dic={}
        for line in fasta_file:
            if line[0]=='>':
                header=line[1:-1]
                header=header[:25]

                fasta_dic[header]=''
            else:
                fasta_dic[header]+=line[:-1]
        fasta_file.close()
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
        fig_with_all_edges=plt.figure()
        fig_with_all_edges.set_figheight(10)
        fig_with_all_edges.set_figwidth(15)
        pos=nx.circular_layout(G)
        weights=[wt for u,v,wt in G.edges(alignments)]
        nx.draw_networkx(G,pos,width=weights)
        labels=nx.get_edge_attributes(G,"weight")
        nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
        plt.axis("on")

        #Create folder
        subprocess.run(["mkdir",tail])

        #Saving all files
        nx.write_gexf(G,f'{tail}/{tail}_with_all_edges.gexf')
        print("Saving the gefx file to {}/{}_with_all_edges.gexf...".format(tail,tail))
        plt.savefig(f"{tail}/{tail}_with_all_edges.png")
        print("Saving the gefx file to {}/{}_with_all_edges.png...".format(tail,tail))

        #normalize the weight
        alignmax,alignmin=max(alignments),min(alignments)
        for i, val in enumerate(alignments):
            alignments[i]=(val-alignmin)/(alignmax-alignmin)

        #Removing all the edges under the threashold, here fixed at 75%
        compt2=0
        for i in range(len(fasta_list_values)):
            for j in range(i+1,len(fasta_list_values)):
                if(alignments[compt2]<0.25):
                    G.remove_edge(fasta_list_keys[i],fasta_list_keys[j])
                compt2+=1




        #Saving all files
        nx.write_gexf(G,f'{tail}/{tail}.gexf')
        print("Saving the gefx file to {}/{}.gexf...".format(tail,tail))
        #plt.savefig(f"{tail}/{tail}.png",dpi=(300))
        #print("Saving the gefx file to {}/{}.png...".format(tail,tail))

        # Button to show the graphic
        self.bouton1 = QPushButton("View the graphic")
        self.bouton1.clicked.connect(self.showing_graph)

        # Button to load Genephi with the file
        self.bouton2 = QPushButton("Load gephi with the file")
        self.bouton2.clicked.connect(self.open_genephi)
 
        # To manage the shape of the window
        layout = QHBoxLayout()
        layout.addWidget(self.bouton1)
        self.setLayout(layout)
        layout.addWidget(self.bouton1)
        layout.addWidget(self.bouton2)
        self.setLayout(layout)
        self.setWindowTitle("Viewing the results of the alignment")

    def showing_graph(self):
        plt.show()
        print("Opening the graphic")

    
    def open_genephi(self): 
        filename=self.fname[0]
        file=filename[:-6]
        home,tail =os.path.split(file)
        gephi_filename="{}/{}_with_all_edges.gexf".format(tail,tail)
        home_path=Path.home()
        path_to_gephi="{}/gephi-0.9.7/bin/gephi".format(home_path)
        result=subprocess.run([path_to_gephi, "-o",gephi_filename]) 
        print("Opening gephi")   




if __name__=="__main__":
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    app = QApplication.instance() 
    if not app: 
        app = QApplication(sys.argv)
    fen = QWidget()
    fen = Fenetre()
    fen.resize(500,250)
    fen.move(300,50)
    fen.show()
    sys.exit(app.exec_())