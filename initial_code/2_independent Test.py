# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 19:54:00 2023

@author: HUAWEI
"""
import numpy as np
import pandas as pd
from bunch import Bunch
#need to run the code in file 1 first
#gene sets at different stages, see Table1 for more details
#choose the gene set/sample you want to use
#RI
sample = ['c_51']#sample for RI model independent Test
gene = ['IGHJ6','IGKC','IGHG4','IGHG3','IGHG1','ANPEP','IGLV3-1']#RI Differential expression
gene = ['GLUL', 'IGKC', 'IGKV4-1', 'FN1', 'JCHAIN', 'TGFBI', 'ACTB', 'SERPINE1', 'CCL19', 'PTGDS',
'VIM','UBC','IGHG4','IGHG2','IGHA1','IGHG1','IGHG3','IGHM','IGHJ6','MT2A','MT1E',
'FTL','IGLV3-1','IGLC1','ANPEP']#RI Chi-square test signatures
gene = ['IGKC','IGHG2','ACTB','IGLV3-1','FTL','IGHA1','IGLC1','JCHAIN','IGHM']#RI final marker

#NRI
sample = ['a_15']#sample for NRI model independent Test
gene = ["IGHG2","IGHG3","IGKC","IGLC2","IGHG1","CD163","IGHGP","IGLC3","IGHG4","JCHAIN","IGHA1","IGLC1","LUM"]#NRI Differential expression
gene = ['APOC1','APOE','B2M','C1QB','CD163','FTL','IGHA1','IGHG1','IGHG2','IGHG3','IGHG4',
 'IGHGP','IGHM','IGKC','IGLC1','IGLC2','IGLC3','JCHAIN','LUM','PSAP','SSR4','TMSB4X']#NRI Chi-square test signatures
gene = ["IGHG3","IGKC","IGLC3","IGHG1","IGHA1","IGHG2"]#NRI final marker

def transformation(Y):#把male,female转化成0,1,数据清洗
    l=[]
    n = 0
    for y in Y:
        if y =='TLS':
            n = n+1
            l.append(1)
        else:
            l.append(0)
    print(n)
    return np.array(l)
def count(Y):
    n = 0
    for y in Y:
        if y ==1:
            n = n+1
    print(n)

TLS_dataset = Bunch()
X = []
y = []

#Independent test
#load the training and validation data
#-------access the annotation from the anootation folder-------
path = 'C:/Users/HUAWEI/Desktop/annotation/c_51.csv'
data = pd.read_csv(path)
c = data.shape[0]
#-------access the preprocessed expression matrix from the matrix folder-------
path = "C:/Users/HUAWEI/Desktop/matrix/c_51.h5 results.csv"
X_Itest = pd.read_csv(path)
X_Itest = X_Itest.loc[gene]
l_gene = X_Itest.index
l_bar = X_Itest.columns
a = X_Itest.shape[0]
b = X_Itest.shape[1]
X_Itest = pd.DataFrame(X_Itest.values.T, columns = l_gene, index = l_bar)
D1 = dict()
key = []
for l in range(c):
    key1 = data['Barcode'][l]
    value1 = data['TLS_2_cat'][l]
    key1 = str(key1).replace('-1','')
    key.append(key1)
    D1[key1] = value1
L1 = list()
for j in l_bar:
    j = str(j).replace('.1','')
    if j in key:
        tls = D1[j]
    else:
        tls = 'NO_TLS'
    L1.append(tls)
y_Itest = transformation(L1)
y_Itest = np.array(y_Itest)
y_Itest = y_Itest.ravel()
count(y_Itest)
#perform the Idependent Test
y_pred2= SVC_.predict(X_Itest)
y_predprob2= SVC_.predict_proba(X_Itest)[:,1]
print("Itest Accuracy : %.4g" % metrics.accuracy_score(y_Itest, y_pred2))
