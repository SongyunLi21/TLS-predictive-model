# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 20:07:49 2023

@author: HUAWEI
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn import metrics
#need to run the code in file 1,2 first
#ROC and PR plot:Figure3/Figure4/FigureS1
#--------ROC--------
#train
pred_y0 = pd.DataFrame(y_pred0)
pred_y0=round(pred_y0,0).astype(int)
fpr0,tpr0,threshold0 = roc_curve(y_train, y_predprob0) ###计算真正率和假正率
roc_auc0 = auc(fpr0,tpr0) ###计算auc的值
#validation
pred_y1 = pd.DataFrame(y_pred1)
pred_y1=round(pred_y1,0).astype(int)
fpr1,tpr1,threshold1 = roc_curve(y_test, y_predprob1) 
roc_auc1 = auc(fpr1,tpr1) 
#Independent test
pred_y2= pd.DataFrame(y_pred2)
pred_y2=round(pred_y2,0).astype(int)
fpr2,tpr2,threshold2 = roc_curve(y_Itest, y_predprob2)
roc_auc2 = auc(fpr2,tpr2) 
#print the result
print("AUROC Score (Train): %f" % metrics.roc_auc_score(y_train, y_predprob0))
print("AUROC Score (Test): %f" % metrics.roc_auc_score(y_test, y_predprob1))
print("AUROC Score (Independent Test): %f" % metrics.roc_auc_score(y_Itest, y_predprob2))
#plot ROC
plt.figure()
plt.rc('font',size=20)
lw = 2
plt.figure(figsize=(10,10))
plt.plot(fpr0, tpr0, color='r',
         lw=lw, label='(ROC={:.4f})(train)'.format(roc_auc0))
plt.plot(fpr1, tpr1, color='g',
         lw=lw, label='(ROC={:.4f})(validation)'.format(roc_auc1))
plt.plot(fpr2, tpr2, color='y',
       lw=lw, label='ROC curve(independent test) (area = %0.2f)' % roc_auc2)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([-0.01, 1.0])
plt.ylim([0.0, 1.01])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic ')
plt.legend(loc="lower right")
#save the plot
sv_path='C:/Users/HUAWEI/Desktop/plots2.0'
plt.savefig(f'{sv_path}/'+'FigureS1A'+'.eps',dpi=600,bbox_inches='tight')
plt.show()
#-----PRC-----
from sklearn.metrics import precision_recall_curve,average_precision_score
#train
precision0, recall0,_ = precision_recall_curve(y_train,y_predprob0)
PRC0 = average_precision_score(y_train,y_predprob0)
area0 = auc(recall0, precision0)
#validation
precision1, recall1, _ = precision_recall_curve(y_test,y_predprob1)
PRC1 = average_precision_score(y_test,y_predprob1)
area1 = auc(recall1, precision1)
#Independent Test
precision2, recall2, _ = precision_recall_curve(y_Itest,y_predprob2)
PRC2 = average_precision_score(y_Itest,y_predprob2)
area2 = auc(recall2, precision2)
#print the result
print("AUPRC Score (Train): %f" % area0)
print("AUPRC Score (Validation): %f" % area1)
print("AUPRC Score (Independent Test): %f" % area2)
#plot PRC
plt.figure(figsize=(10,10))
plt.step(recall0, precision0, color='r', label=' (PRC={:.4f})(train)'.format(area0))
plt.step(recall1, precision1, color='g', label=' (PRC={:.4f})(validation)'.format(area1))
plt.step(recall2, precision2, color='y', label=' (PRC={:.4f})(Independent test)'.format(area2)) 
plt.plot([0, 1], [1, 0], color='m', linestyle='--')
plt.title('PR curves')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision Recall Curves')
plt.legend(loc='lower right')
#save the plot
sv_path='C:/Users/HUAWEI/Desktop/plots2.0'
plt.savefig(f'{sv_path}/'+'FigureS1B'+'.eps',dpi=600,bbox_inches='tight')
plt.show()

#------------spatial visulization-------------------------
sample = ['c_3','c_4','c_36','c_51']#RI
sample = ['a_3','a_15','b_1','b_18']#NRI
L_data2 = []
L_bartls = []
for s in range(len(sample)):
    #-------access the preprocessed expression matrix from the matrix folder-------
    path = "C:/Users/HUAWEI/Desktop/matrix/"+sample[s]+".h5 results.csv"
    data1 = pd.read_csv(path)
    data1 = data1.loc[gene]#when using the original model.Delete this line
    l_gene = data1.index
    l_barcodes1 = data1.columns
    l_bar = []
    for i in l_barcodes1:
        i=i.replace('.1','-1')
        l_bar.append(i)
    a = data1.shape[0]
    b = data1.shape[1]
    data1 = pd.DataFrame(data1.values.T, columns = l_gene, index = l_bar)
    data1 = SVC_.predict(data1)
    data1 = pd.DataFrame(data1)
    data1 = pd.DataFrame(data1.values, index = l_bar)
    #----------------access from the postition folder--------------------
    path2 = 'C:/Users/HUAWEI/Desktop/position/'+sample[s]+'_positions_list.csv'
    data2 = pd.read_csv(path2)
    data2 = data2[['barcodes','X','Y']]
    data2 = pd.DataFrame(data2.values.T,index = ['barcodes','X','Y'])
    l_barcodes2 = data2.loc['barcodes'].tolist()
    data2 = pd.DataFrame(data2.values,columns = data2.loc['barcodes'])
    data2.index = ['barcodes','X','Y']
    #connect position and annotation
    D_tls = {}
    for key in l_bar:
        D_tls[key] = data1.loc[key]
    l_bartls = []
    for bar in l_barcodes2:
        if bar in l_bar:
            l_bartls.append(int(D_tls[bar]))
        else:
            l_bartls.append(-1)
    L_data2.append(data2)
    L_bartls.append(l_bartls)
#spatial plot
import matplotlib as mpl
cmp = mpl.colors.ListedColormap(['lavender','darkblue','gold'])
norm = mpl.colors.BoundaryNorm([-1,0,1,2], cmp.N)
fig, ax = plt.subplots(1, 4, figsize=(24,5))
for i in range(len(sample)):
    x = 0
    y = i
    name = sample[i]
    ax[y].set_title(name, fontsize=30)
    ax[y].get_xaxis().set_visible(False)
    ax[y].get_yaxis().set_visible(False)
    ax[y].scatter(np.array(L_data2[i].loc['X'].tolist()),np.array(L_data2[i].loc['Y'].tolist()), c=L_bartls[i], edgecolors='none', cmap=cmp, norm=norm,alpha=0.7)
#save the plot
sv_path='C:/Users/HUAWEI/Desktop/plots2.0'
plt.savefig(f'{sv_path}/'+'Figure4A'+'.eps',dpi=600,bbox_inches='tight')#保存文件在指定文件夹下很方便
plt.show()



