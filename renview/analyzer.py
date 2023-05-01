#Dondi 2023

import numpy as np

def NodeAnalyzer(nodenum,code,path): #1 R c:\data
    edgein=[]
    edgeinvalues=[]
    edgeout=[]
    edgeoutvalues=[]
    edgestokeep=[]
    nodestokeep=[]
    filename=path+str(nodenum)+".txt"
    nodefile = open(filename)
    Lines = nodefile.readlines()
    for line in Lines:
        if 'edge[' in line:
            temp=line
            e=temp.find('mol')-1
            s=temp.find('%')+1
            tempvalue=float(temp[s:e])
        if '->' in line:
            if '"'+code+str(nodenum)+'"-' in line: #compound indicated in nodenum is a reactant
                edgein.append(temp)
                edgeinvalues.append(tempvalue)
            else:  #compound indicated in nodenum is a reactant
                edgeout.append(temp)
                edgeoutvalues.append(tempvalue)
    totin=np.sum(edgeinvalues)            
    avgin=np.average(edgeinvalues)
    stddevin=np.std(edgeinvalues)
    print("Total input speed=",totin)
    print(" average=",avgin)    
    print(" std deviation=",stddevin)
    if stddevin<0.2*avgin:
        print(" --- Almost no variation")
    print()
    totout=np.sum(edgeoutvalues)
    avgout=np.average(edgeoutvalues)
    stddevout=np.std(edgeoutvalues)
    print("Total output speed=",totout)
    print(" average=",avgout)    
    print(" std deviation=",stddevout)
    if stddevout<0.2*avgout:
        print("Almost no variation")    

NodeAnalyzer(1,'R','./results/example_caz/Species/')            
    
