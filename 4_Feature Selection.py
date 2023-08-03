# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 20:19:52 2023

@author: HUAWEI
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sns
#the code in this file is optional
#Files 1-3 contain all the process for obtaining the final model and most of the figures
#The code in this file should only be used if you want to verify the feature selection process. 
#--------------------------------------------------------------------------------------
#gene sets at different stages, see Table1 for more details
#choose the gene set/sample you want to use
#RI
gene = ['IGHJ6','IGKC','IGHG4','IGHG3','IGHG1','ANPEP','IGLV3-1']#RI Differential expression
gene = ['GLUL', 'IGKC', 'IGKV4-1', 'FN1', 'JCHAIN', 'TGFBI', 'ACTB', 'SERPINE1', 'CCL19', 'PTGDS',
'VIM','UBC','IGHG4','IGHG2','IGHA1','IGHG1','IGHG3','IGHM','IGHJ6','MT2A','MT1E',
'FTL','IGLV3-1','IGLC1','ANPEP']#RI Chi-square test signatures
gene = ['IGKC','IGHG2','ACTB','IGLV3-1','FTL','IGHA1','IGLC1','JCHAIN','IGHM']#RI final marker

#NRI
gene = ["IGHG2","IGHG3","IGKC","IGLC2","IGHG1","CD163","IGHGP","IGLC3","IGHG4","JCHAIN","IGHA1","IGLC1","LUM"]#NRI Differential expression
gene = ['APOC1','APOE','B2M','C1QB','CD163','FTL','IGHA1','IGHG1','IGHG2','IGHG3','IGHG4',
 'IGHGP','IGHM','IGKC','IGLC1','IGLC2','IGLC3','JCHAIN','LUM','PSAP','SSR4','TMSB4X']#NRI Chi-square test signatures
gene = ["IGHG3","IGKC","IGLC3","IGHG1","IGHA1","IGHG2"]#NRI final marker
#------------------model selection--------------------
#this code should be run after the data in file 1 is loaded and before the model construction in file 1
from sklearn.linear_model import LogisticRegression as LR
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score,StratifiedKFold
from sklearn.tree import DecisionTreeClassifier as DTC
LR_ = LR()
GB_ = GaussianNB()
SVC_ = SVC()
DTC_ = DTC()

method = [LR_,GB_,SVC_,DTC_]
name = ['LogisticRegression','GaussianNB','SVM Classifier','DecisionTreeClassifier']
for j in range(len(method)):
    m = method[j]
    strKFold = StratifiedKFold(n_splits=3,shuffle=False)
    scores = cross_val_score(m,X,y,cv=strKFold)
    print("---------{}--------".format(name[j]))
    print("leave-one-out cross validation scores of:{}".format(scores))
    print("Mean score of leave-one-out cross validation:{:.2f}".format(scores.mean()))
    
#------------------Chi Square Test---------------------------------
#this code should be run after the data in file 1 is loaded and before the model construction in file 1
from sklearn.feature_selection import chi2
from sklearn.feature_selection import SelectKBest
X_fsvar = VarianceThreshold(np.median(TLS_dataset.data.var().values)).fit_transform(TLS_dataset.data)
score = []
for i in range(10,30,2):#the range can be change, see Figure S2
    X_fschi = SelectKBest(chi2, k=i).fit_transform(X_fsvar, y)
    once = cross_val_score(SVC_,X_fschi,TLS_dataset.target,cv=5).mean()
    print(once)
    score.append(once)
plt.plot(range(10,30,2),score,size = 10)
plt.title('chi square test')
plt.xlabel('number of genes')
plt.ylabel('accuracy')
plt.title('chi square test',fronsize = 20)
sv_path='C:/Users/HUAWEI/Desktop/plots2.0'
plt.savefig(f'{sv_path}/'+'FigureS2D'+'.eps',dpi=600,bbox_inches='tight')
plt.show()
X_fschi = SelectKBest(chi2, k=24).fit_transform(TLS_dataset.data, TLS_dataset.target)   

#------------------permutation Importance---------------------------
#this code should be run after the model construction in file 1 is finished
from sklearn.inspection import permutation_importance
result = permutation_importance(SVC_, X_test, y_test)
featImp = pd.DataFrame()
featImp['feat'] = l_gene
featImp['importance'] = result.importances_mean
featImp = featImp.sort_values('importance',ascending = False)
plt.figure(figsize=[15,10])
plt.title("Permutation importance",fronsize = 30)
sns.barplot(x = 'importance', y = 'feat',data = featImp[:len(gene)],orient='h')
sv_path='C:/Users/HUAWEI/Desktop/plots2.0'
plt.savefig(f'{sv_path}/'+'Figure2'+'.eps',dpi=600,bbox_inches='tight')
plt.show()
