# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 08:51:08 2018

@author: mlee
"""

#%% dictionary of KEGG entrezID2pathwayID and pathwayID2entrezID
fKEGG="D:\miRNADB\work_python\miRNetwork\KEGG\pid2gid";

dKEGG_geneKEY={};
f=open(fKEGG,'r');
for line in f:
    line=line.split();
    key=line[1][4:];
    val=line[0];
    
    vals=dKEGG_geneKEY.get(key);
    if bool(vals)==False:        
        dKEGG_geneKEY[key]=[val];
    else:        
        vals.append(val);
        dKEGG_geneKEY[key]=list(set(vals));
    
f.close()   
len(dKEGG_geneKEY.keys())


dKEGG_pathwayKEY={};
f=open(fKEGG,'r');
for line in f:
    line=line.split();
    key=line[0];
    val=line[1][4:];
    
    vals=dKEGG_pathwayKEY.get(key);
    if bool(vals)==False:        
        dKEGG_pathwayKEY[key]=[val];
    else:        
        vals.append(val);
        dKEGG_pathwayKEY[key]=list(set(vals));
    
f.close()   
len(dKEGG_pathwayKEY.keys())


#%% dictionary of miRTarBase miR2symbol and symbol2miR
fmiRTarBase="D:\miRNADB\work_python\miRNetwork\miRTarBase\hsa_MTI_sum.txt";

dmiRTarBase_miRKEY={};
f=open(fmiRTarBase,'r');
for line in f:
    line=line.split();
    key=line[0];
    val=line[1];
    
    vals=dmiRTarBase_miRKEY.get(key);
    if bool(vals)==False:
        dmiRTarBase_miRKEY[key]=[val];
    else:        
        vals.append(val);
        dmiRTarBase_miRKEY[key]=list(set(vals));
    
f.close() 
del dmiRTarBase_miRKEY['miRNA'];
len(dmiRTarBase_miRKEY.keys())   


dmiRTarBase_geneKEY={};
f=open(fmiRTarBase,'r');
for line in f:
    line=line.split();
    key=line[1];
    val=line[0];
    
    vals=dmiRTarBase_geneKEY.get(key);
    if bool(vals)==False:
        dmiRTarBase_geneKEY[key]=[val];
    else:        
        vals.append(val);
        dmiRTarBase_geneKEY[key]=list(set(vals));
    
f.close() 
del dmiRTarBase_geneKEY['target'];
len(dmiRTarBase_geneKEY.keys()) 


#%% dictionary of NCBI symbol2entrezID and entrezID2symbol    
fNCBI="D:\miRNADB\work_python\miRNetwork\\NCBI\hsa_gene_info";

dNCBI_symbolKEY={};
f=open(fNCBI,'r');
for line in f:
    line=line.split();
    key=line[2];
    val=line[1];
    
    vals=dNCBI_symbolKEY.get(key);
    if bool(vals)==False:
        dNCBI_symbolKEY[key]=[val];
    else:        
        vals.append(val);
        dNCBI_symbolKEY[key]=list(set(vals));
    
f.close() 
del dNCBI_symbolKEY['Symbol'];
len(dNCBI_symbolKEY.keys())   


dNCBI_entrezIDKEY={};
f=open(fNCBI,'r');
for line in f:
    line=line.split();
    key=line[1];
    val=line[2];
    
    vals=dNCBI_entrezIDKEY.get(key);
    if bool(vals)==False:
        dNCBI_entrezIDKEY[key]=[val];
    else:        
        vals.append(val);
        dNCBI_entrezIDKEY[key]=list(set(vals));
    
f.close() 
del dNCBI_entrezIDKEY['GeneID'];
len(dNCBI_entrezIDKEY.keys())   


#%% dictionary of miR2pathwayID and pathwayID2miR
miR=list(dmiRTarBase_miRKEY.keys());
miR.sort();

miR2KEGG={};
for x in miR:
    symbol=dmiRTarBase_miRKEY.get(x);
    
    pathways=[];
    for i in range(0,len(symbol)):
        entrezID=dNCBI_symbolKEY.get(symbol[i]);
        if bool(entrezID):
            pathway=dKEGG_geneKEY.get(entrezID[0]);
            if bool(pathway):
                pathways=pathways+pathway;        
    miR2KEGG[x]=list(set(pathways));


KEGG=list(dKEGG_pathwayKEY.keys());
KEGG.sort();
KEGG2miR={};
for x in KEGG:
    entrezID=dKEGG_pathwayKEY.get(x);
    
    miRs=[];
    for i in range(0,len(entrezID)):
        symbol=dNCBI_entrezIDKEY.get(entrezID[i]);
        if bool(symbol):
            miR=dmiRTarBase_geneKEY.get(symbol[0]);
            if bool(miR):
                miRs=miRs+miR;        
    KEGG2miR[x]=list(set(miRs));


#%% fisher's exact test
f=open("DEmiR.txt",'r');
DEmiR=f.read().splitlines();
f.close();

miR=list(miR2KEGG.keys());
nonDEmiR=list(set(miR) - set(DEmiR));

import numpy as np
import scipy.stats as stats

statFisher=np.zeros((len(KEGG),2));
i=0;
for pathway in KEGG:
    miR_in_pathway=set(KEGG2miR.get(pathway));
    miR_NOTin_pathway=set(miR) - set(miR_in_pathway);
    
    a=set(DEmiR) & set(miR_in_pathway);
    b=set(DEmiR) & set(miR_NOTin_pathway);
    c=set(nonDEmiR) & set(miR_in_pathway);
    d=set(nonDEmiR) & set(miR_NOTin_pathway);
    obs = np.array([[len(a),len(b)], [len(c),len(d)]]);
    odd, pval = stats.fisher_exact(obs);
    statFisher[i,0]=odd;
    statFisher[i,1]=pval;
    i=i+1;
    

import pandas as pd

data1=pd.DataFrame(KEGG)
data2=pd.DataFrame(statFisher)
data=pd.concat([data1,data2], axis=1)
data.columns=['KEGG','Odds ratio','P-value'];
data.to_csv('stats_FisherExactTest.txt', header=True, index=True, sep='\t')
















    
        


