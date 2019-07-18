#random forrest
import matplotlib
matplotlib.use('Agg')
from sklearn.utils import shuffle
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from sklearn.decomposition import PCA
from sklearn.utils import shuffle
from sklearn.preprocessing import normalize
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from collections import defaultdict
import operator
import matplotlib.pyplot as plt 
import random
from joblib import Parallel, delayed
import multiprocessing
from sklearn import datasets

from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

from sklearn.ensemble import BaggingClassifier


import itertools
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB 
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import confusion_matrix
from sklearn.feature_selection import SelectFromModel

from sklearn.metrics import precision_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.model_selection import train_test_split
from random import randint
from mlxtend.classifier import StackingClassifier
import sys
import warnings
warnings.filterwarnings("ignore")
from sklearn.model_selection import cross_val_score, train_test_split


def random_forrest(cancer_type):
    feat_type=[]
    data_dir = "/home/ubuntu/cancer/"
    data_file = data_dir + cancer_type + "_matrix.csv"
    output_file = data_dir + cancer_type + "_output.txt"
    graph_file = data_dir + cancer_type + "_graph.png"
    file = open(output_file, "w")
    df = pd.read_csv(data_file)
    df_normal=df[df['label']==0]
    df_tumor=df[df['label']==1]
    unique, counts = np.unique(df['label'], return_counts=True)
    f=dict(zip(unique, counts))
    features=[]
    rf = RandomForestRegressor()
    scores = defaultdict(list)
    #class imbalance
    for z in range(500):
        q=random.randint(0,(df_tumor.shape[0]-df_normal.shape[0]))
        frames=[df_normal,df_tumor[q:q+df_normal.shape[0]]]
        result=pd.concat(frames)
        result = shuffle(result)
        file_ids=result.pop('file_id')
        Y=result.pop('label').values
        names=result.columns
        a=normalize(result)
        dafr=pd.DataFrame(a)
        X = a
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=0)
        r = rf.fit(X_train, y_train)
        acc = r2_score(y_test, rf.predict(X_test))
        for i in range(X.shape[1]):
            X_t = X_test.copy()
            np.random.shuffle(X_t[:, i])
            shuff_acc = r2_score(y_test, rf.predict(X_t))
            scores[names[i]].append((acc-shuff_acc)/acc)
        feat=sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), names), reverse=True)
        features.append(feat)
    #averaging important features
    avg={}
    for m in range(25):
        c=[]
        for j in range(10):
            for l in range(50):
                if features[j][l][1]==features[0][m][1]:
                    break
            c.append(features[j][l][0])
        avg[features[0][m][1]]=np.mean(c)
        sorted_feat = sorted(avg.items(), key=operator.itemgetter(1), reverse=True)
    feat_type.append(sorted_feat)
    feat_value=[]
    name=[]
    file.write('Feature Importances'+'\n')
    #20 most important features
    for g in range(20):
        feat_value.append(sorted_feat[g][1])
        name.append(sorted_feat[g][0])
        file.write(sorted_feat[g][0]+'\n')
    #important features plot
    plt.figure()
    plt.title('Feature Importances')
    feat_importances = pd.Series(feat_value, index=name)
    feat_importances.nlargest(10).plot(kind='barh')
    plt.savefig(graph_file)
    file.close()   


def data_ensemble(cancer_type,feat):
	data_dir = "/home/ubuntu/cancer/"
	data_file = data_dir + cancer_type + "_matrix.csv"
	features = data_dir + cancer_type + "_output.txt"
	output_file = data_dir + cancer_type + "_accuracy.txt"
	file = open(features, "r")
	o_file = open(output_file, "w")
	line = file.readline()
	line = file.readline()
	df = pd.read_csv(data_file)
	df = shuffle(df)
	file_ids=df.pop('file_id')
	y = df.pop('label').values
	dataf=df.pop(line[:-1])
	#dataframe consisting of only important features
	for x in range(feat):
		line = file.readline()
		dataf=np.column_stack((dataf,df.pop(line[:-1])))
	X=normalize(dataf)
	X=scale(X)
	pca=PCA()
	pca.fit(X)
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42)
        #multiple classifiers
	clf1 = RandomForestClassifier(random_state=1,n_estimators=100)
	clf2 = GradientBoostingClassifier(n_estimators=1200,subsample=0.5,random_state=3)
	clf3 = SVC(gamma='auto')
	clf4 = KNeighborsClassifier(n_neighbors=1)
	clf5 = DecisionTreeClassifier(random_state=0)
	lr = LogisticRegression(solver='lbfgs')
	#stacking for data ensemble
	sclf = StackingClassifier(classifiers=[clf1, clf2, clf3, clf4, clf5], meta_classifier=lr)
	clf1.fit(X_train,y_train)
	clf2.fit(X_train,y_train)
	clf3.fit(X_train,y_train)
	clf4.fit(X_train,y_train)
	clf5.fit(X_train,y_train)
	sclf.fit(X_train,y_train)
	y_test_predict=sclf.predict(X_test)
	precision = precision_score(y_test, y_test_predict)
	accuracy = accuracy_score(y_test, y_test_predict)
	f1 = f1_score(y_test, y_test_predict)
	recall = recall_score(y_test, y_test_predict)
	scores = [precision,accuracy,f1,recall]
	label = ['RF', 'GBDT', 'SVM','KNN','DT','Stacking']
	clf_list = [clf1, clf2, clf3, clf4, clf5, sclf]
	#score calculation
	for clf, label in zip(clf_list, label):
		y_test_predict = clf.predict(X_test)
		tn, fp, fn, tp = confusion_matrix(y_test, y_test_predict).ravel()
		specificity = tn / (tn+fp)
		recall = tp / (tp+fn)
		precision = tp / (tp+fp)
		accuracy = (tp + tn) / (tp+tn+fp+fn)
		f1 = 2*tp / (2*tp+fp+fn)
		o_file.write("\nAccuracy: %.2f [%s] \nPrecision: %.2f [%s] \nRecall: %.2f [%s] \nF1 score: %.2f [%s] \nSpecificity: %.2f [%s]\n" %(accuracy,label,precision, label, recall, label, f1, label, specificity, label))

if __name__ == '__main__':
    cancer = {"Lung" : 7, "Breast" : 5, "Kidney" : 7, "Stomach" : 10, "Head and Neck" : 7, "Prostate" : 18, "Thyroid" : 7, "Liver" : 10}
    a = datetime.datetime.now().replace(microsecond=0)
    lst=["Breast","Lung","Liver","Kidney","Prostate","Head and Neck","Stomach","Thyroid"]
    num_cores = multiprocessing.cpu_count()
    rf = Parallel(n_jobs=num_cores, verbose=50)(delayed(random_forrest)(cancer_type)for cancer_type in lst)
    data_en = Parallel(n_jobs=num_cores, verbose=50)(delayed(data_ensemble)(cancer_type,cancer[cancer_type])for cancer_type in lst)
    b = datetime.datetime.now().replace(microsecond=0)
    #code run time
    print("Code computation time(hrs:min:sec): ")
    print(b-a)
