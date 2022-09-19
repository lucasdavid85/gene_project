# importations à faire pour la réalisation d'une interface graphique
import sys

from PyQt5.QtWidgets import QApplication, QWidget,QPushButton,QHBoxLayout,QFileDialog

#Library importation
import os.path
import networkx as nx 
from Bio import pairwise2
import matplotlib.pyplot as plt
import sys
import subprocess
import os.path
from pathlib import Path




class Fenetre(QWidget):

    def __init__(self):
        QWidget.__init__(self)
#
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
        pos=nx.circular_layout(G)
        labels=nx.get_edge_attributes(G,'weight')
        nx.draw(G,pos,with_labels=True,node_size=15, font_size = 10,node_color="blue",width=0.5)
        esmall = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] <= 0.5]
        nx.draw_networkx_edges(G, pos, edgelist=esmall, width=10, alpha=0.5, edge_color="b")
        nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
        plt.axis("on")

        #Create folder
        subprocess.run(["mkdir",tail])

        #Saving all files
        nx.write_gexf(G,f'{tail}/{tail}.gexf')
        print("Saving the gefx file to {}/{}.gexf...".format(tail,tail))
        plt.savefig(f"{tail}/{tail}.png")
        print("Saving the gefx file to {}/{}.png...".format(tail,tail))



        # creation du premier bouton
        self.bouton1 = QPushButton("View the graphic")
        self.bouton1.clicked.connect(self.showing_graph)


        # creation du deuxieme bouton
        self.bouton2 = QPushButton("Load gephi with the file")
        self.bouton2.clicked.connect(self.open_genephi)
 
        #self.bouton3 = QPushButton("Entrer un nouveau fichier fasta")
        #self.bouton3.clicked.connect(self.upload_fasta)


        # creation du gestionnaire de mise en forme de type QHBoxLayout
        layout = QHBoxLayout()
        layout.addWidget(self.bouton1)
        self.setLayout(layout)


        # ajout du premier bouton au gestionnaire de mise en forme
        layout.addWidget(self.bouton1)
        # ajout du deuxieme bouton au gestionnaire de mise en forme
        layout.addWidget(self.bouton2)
        #layout.addWidget(self.bouton3)
        # on fixe le gestionnaire de mise en forme de la fenetre

        self.setLayout(layout)

        self.setWindowTitle("Ma fenetre")



    #def upload_fasta(self): 
    #    dialog = QFileDialog()
    #    self.fname = dialog.getOpenFileName(None, "Import FASTA", "", "FASTA data files (*.fasta)")


    def showing_graph(self):
        plt.show()
        print("Opening the graphic")

       

    def open_genephi(self): 
        filename=self.fname[0]
        file=filename[:-6]
        home,tail =os.path.split(file)
        gephi_filename="{}/{}.gexf".format(tail,tail)
        home_path=Path.home()
        path_to_gephi="{}/gephi-0.9.7/bin/gephi".format(home_path)
        result=subprocess.run([path_to_gephi, "-o",gephi_filename]) 
        print("Opening gephi")   







         
app = QApplication.instance() 

if not app: # sinon on crée une instance de QApplication
    app = QApplication(sys.argv)


fen = QWidget()
fen = Fenetre()
#tail=fen.upload_fasta()
#fen.open_genephi(tail)
fen.resize(500,250)
fen.move(300,50)
fen.show()


#w = MainWindow()
#w.show()

sys.exit(app.exec_())