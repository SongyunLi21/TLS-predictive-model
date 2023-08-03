# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 20:17:01 2023

@author: HUAWEI
"""
import numpy as np
import pandas as pd
from bunch import Bunch

#gene sets at different stages, see Table1 for more details
#choose the gene set/sample you want to use
#RI
sample = ['c_3','c_4','c_36']#the samples used for trainning and validation of RI model
gene = ['IGHJ6','IGKC','IGHG4','IGHG3','IGHG1','ANPEP','IGLV3-1']#RI Differential expression
gene = ['GLUL', 'IGKC', 'IGKV4-1', 'FN1', 'JCHAIN', 'TGFBI', 'ACTB', 'SERPINE1', 'CCL19', 'PTGDS',
'VIM','UBC','IGHG4','IGHG2','IGHA1','IGHG1','IGHG3','IGHM','IGHJ6','MT2A','MT1E',
'FTL','IGLV3-1','IGLC1','ANPEP']#RI Chi-square test signatures
gene = ['IGKC','IGHG2','ACTB','IGLV3-1','FTL','IGHA1','IGLC1','JCHAIN','IGHM']#RI final marker

#NRI
sample = ['a_3','b_1','b_18']#the samples used for trainning and validation of NRI model
gene = ["IGHG2","IGHG3","IGKC","IGLC2","IGHG1","CD163","IGHGP","IGLC3","IGHG4","JCHAIN","IGHA1","IGLC1","LUM"]#NRI Differential expression
gene = ['APOC1','APOE','B2M','C1QB','CD163','FTL','IGHA1','IGHG1','IGHG2','IGHG3','IGHG4',
 'IGHGP','IGHM','IGKC','IGLC1','IGLC2','IGLC3','JCHAIN','LUM','PSAP','SSR4','TMSB4X']#NRI Chi-square test signatures
gene = ["IGHG3","IGKC","IGLC3","IGHG1","IGHA1","IGHG2"]#NRI final marker

def transformation(Y):
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

#load the training and validation data
for i in range(len(sample)):
    #-------access the annotation from the anootation folder-------
    path = 'C:/Users/HUAWEI/Desktop/annotation/'+sample[i]+'.csv'
    data = pd.read_csv(path)
    c = data.shape[0]
    #-------access the preprocessed expression matrix from the matrix folder-------
    path = "C:/Users/HUAWEI/Desktop/matrix/"+sample[i]+".h5 results.csv"
    X_data = pd.read_csv(path)
    X_data = X_data.loc[gene]#for original model, delete this line
    l_gene = X_data.index
    l_bar = X_data.columns
    a = X_data.shape[0]
    b = X_data.shape[1]
    X_data = pd.DataFrame(X_data.values.T, columns = l_gene, index = l_bar)
    D1 = dict()
    key = []
    for l in range(c):
        x = '-'+str(i)
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
    y_data = transformation(L1)
    y_data = y_data.ravel()
    X_data = np.array(X_data).tolist()
    y_data = np.array(y_data).tolist()
    X = X + X_data
    y = y + y_data
X = pd.DataFrame(X)
X = np.array(X)
y = np.array(y)
TLS_dataset.data = pd.DataFrame(X,columns = l_gene)
TLS_dataset.target = y

#------------------Model Construction-----------------------
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
X_train, X_test, y_train, y_test = train_test_split(
    TLS_dataset.data, TLS_dataset.target, random_state = 0)
count(y_train)
count(y_test)
SVC_ = SVC(probability=True, class_weight={1:0.8,0:0.2})
#RI final model: class_weight={1:0.8,0:0.2}
#NRI final model: class_weight={1:0.9,0:0.1}
#see TablecS4 for more details about the class weight information
SVC_.fit(X_train,y_train)
from sklearn import metrics
y_pred0= SVC_.predict(X_train)
y_predprob0= SVC_.predict_proba(X_train)[:,1]
print("Training Accuracy : %.4g" % metrics.accuracy_score(y_train, y_pred0))
y_pred1= SVC_.predict(X_test)
y_predprob1= SVC_.predict_proba(X_test)[:,1]
print("Testing Accuracy : %.4g" % metrics.accuracy_score(y_test, y_pred1))