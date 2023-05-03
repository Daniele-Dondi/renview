#Dondi 2023
#C:\Users\utente\AppData\Local\Microsoft\WindowsApps\PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0
#C:\Users\utente\AppData\Local\Microsoft\WindowsApps

import numpy as np
from graphviz import *
from subprocess import check_call

def NodeAnalyzer(nodenum,code,path,percentile,prefix): #1 R c:\data
    edgein=[]
    edgeinvalues=[]
    edgeintokeep=[]
    edgeout=[]
    edgeoutvalues=[]
    edgeouttokeep=[]
    nodestokeep=[]
    filename=path+str(nodenum)+".txt"
    fnameout=path+prefix+str(nodenum)+".txt"
    svgout=path+prefix+str(nodenum)+".svg"
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
                edgeout.append(temp+line)
                edgeoutvalues.append(tempvalue)
            else:  #compound indicated in nodenum is a product
                edgein.append(temp+line)
                edgeinvalues.append(tempvalue)
    totin=np.sum(edgeinvalues)            
    avgin=np.average(edgeinvalues)
    stddevin=np.std(edgeinvalues)
    medianin=np.median(edgeinvalues)
    percentilein=np.percentile(edgeinvalues,percentile)
    print("Node ",nodenum," analysis:\n")
    print("Total number of inputs:",len(edgeinvalues))
    print(" Total input speed=",totin)    
    print(" average=",avgin)    
    print(" std deviation=",stddevin)
    if stddevin<0.2*avgin:
        print(" --- Almost no variation")
    print(" median=",medianin)
    print(" "+str(percentile)+"th percentile=",percentilein)    
    for idx, x in enumerate(edgeinvalues):
        if x>=percentilein:
            edgeintokeep.append(edgein[idx])
    print(" Reduced number of inputs:",len(edgeintokeep))            
    print()
    totout=np.sum(edgeoutvalues)
    avgout=np.average(edgeoutvalues)
    stddevout=np.std(edgeoutvalues)
    medianout=np.median(edgeoutvalues)
    percentileout=np.percentile(edgeoutvalues,percentile)    
    print("Total number of otputs:",len(edgeoutvalues))
    print(" Total output speed=",totout)    
    print(" average=",avgout)    
    print(" std deviation=",stddevout)
    if stddevout<0.2*avgout:
        print(" --- Almost no variation")
    print(" median=",medianout)
    print(" "+str(percentile)+"th percentile=",percentileout)   
    for idx, x in enumerate(edgeoutvalues):
        if x>=percentileout:
            edgeouttokeep.append(edgeout[idx])
    print(" Reduced number of outputs:",len(edgeouttokeep))                        
    #calculates nodes to keep
    for idx, x in enumerate(edgeintokeep):
        x=x.splitlines()[1].split('"')
        if x[1] not in nodestokeep:
            nodestokeep.append(x[1])
    for idx, x in enumerate(edgeouttokeep):
        x=x.splitlines()[1].split('"')
        if x[3] not in nodestokeep:
            nodestokeep.append(x[3])            
    print("\nNodes kept: ",len(nodestokeep)," out of ",len(edgeinvalues)+len(edgeoutvalues))
    nodestokeep.append(code+str(nodenum))
    
    f = open(fnameout, "w+")
    f.write('digraph G {\n')
    f.write('splines = true;\n')
    f.write('graph [bgcolor=lightgray, resolution=64, fontname=Arial, fontcolor=blue, fontsize=36];\n')
    f.write('node [fontsize=12];\n')
    f.write('edge [fontsize=45];\n') #set312,, colorscheme=paired12
    f.write('label = "Reaction Path Analysis";\n')
    f.write('labelloc = "t";\n')
    f.write('center=1;\n')
    f.write('size="10,10";\n')
    f.write('ranksep="0.25 equally";\n')
    f.write('nodesep="0.25 equally";\n')
    f.write('rankdir=LR;\n')
    f.write('bgcolor=white;\n')
    for line in Lines:
        if '"[' in line:
            x=line.split('"')
            if (x[1]) in nodestokeep:
                f.write(line)
    for line in edgeouttokeep:
        f.write(line)
    for line in edgeintokeep:
        f.write(line)        
    f.write('}\n')
    f.write('\n')    
    f.close()
    check_call(['dot', '-Tsvg', fnameout, '-o',svgout])
    
NodeAnalyzer(13,'R','./results/example_caz/Species/',99,'simplyfied-')            
    
