#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:25:19 2019

@author: jun
"""
import os
os.chdir("/home/jun/Desktop/lncRNAs_ASD/Dataset/feature_selected/")

import os,sys
import getopt
import pandas as pd
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix, matthews_corrcoef
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from keras.models import Sequential, Model
from keras.layers import Input, Dense, Activation, Dropout
from keras import regularizers
from keras import optimizers
from sklearn import preprocessing


def build_deep_classical_autoencoder(input_dim, encoding_dim):
    input_fea = Input(shape=(input_dim,))
    encoded = Dense(encoding_dim, activation='relu')(input_fea)
    #encoded = Dense(int(hidden_dim/2), activation='relu')(encoded)
    #encoded = Dense(int(hidden_dim/4), activation='relu')(encoded)
    #decoded = Dense(int(hidden_dim/2), activation='relu')(encoded)
    #decoded = Dense(hidden_dim, activation='relu')(encoded)
    decoded = Dense(input_dim, activation='relu')(encoded)
    autoencoder = Model(input_fea, decoded)
    encoder = Model(input_fea, encoded)
    
    return autoencoder, encoder

def deep_autoencoder(X):
    autoencoder, encoder = build_deep_classical_autoencoder(524, 16) 
    autoencoder.compile(loss='mean_squared_error', optimizer='RMSprop')
    autoencoder.fit(X, X, batch_size = 64, epochs = 550, shuffle=False)
    encoded_X = encoder.predict(X, verbose=0)
    print(encoded_X.shape)
    return encoded_X

def para_tuning(model,para,X,Y):
    score = make_scorer(strength, greater_is_better=True)
    grid_obj = GridSearchCV(model, para, cv=5, scoring = score)
    grid_obj = grid_obj.fit(X, Y)
    para_best = grid_obj.best_estimator_
    return(para_best)

def strength(Ytrue,Ypred):
    confusion = confusion_matrix(Ytrue, Ypred)
    TP = confusion[1, 1]
    TN = confusion[0, 0]
    FP = confusion[0, 1]
    FN = confusion[1, 0]
    sepcificity = TN / float( TN + FP)
    sensitivity = TP / float(FN + TP)
    strengthv= (sepcificity+sensitivity)/2
    return(strengthv)

def cv(model, X_train, Y_train, random_seed, fileout):
    seed = int(random_seed)
    cv_outcomes_acc = []
    cv_outcomes_sep = []
    cv_outcomes_sen = []
    cv_outcomes_str = []
    cv_outcomes_mcc = []
    cv_outcomes_f1 = []
    cv_outcomes_auc = []   
    cv_preds = []
    cv_probs = []
    kf = model_selection.KFold(n_splits=10, shuffle = True, random_state=seed)
    for train_index, test_index in kf.split(X_train):
        X_cv_train, X_cv_test = X_train[train_index,:], X_train[test_index,:]
        Y_cv_train, Y_cv_test = Y_train[train_index], Y_train[test_index]
        model.fit(X_cv_train,Y_cv_train)
        predictions = model.predict(X_cv_test)
        pred_train_prob = model.predict_proba(X_cv_test)[:, 1]
        accuracy = accuracy_score(Y_cv_test,predictions)
        confusion = confusion_matrix(Y_cv_test, predictions)
        TP = confusion[1, 1]
        TN = confusion[0, 0]
        FP = confusion[0, 1]
        FN = confusion[1, 0]
        sepcificity = TN / float( TN + FP)
        sensitivity = TP / float(FN + TP)
        strength = (sepcificity + sensitivity)/2
        mcc = matthews_corrcoef(Y_cv_test, predictions)
        f1 = f1_score(Y_cv_test, predictions)
        fpr, tpr, thresholds = roc_curve(Y_cv_test, pred_train_prob)
        cv_preds = cv_preds + Y_cv_test.tolist()
        cv_probs = cv_probs + pred_train_prob.tolist()
        aucvalue = auc(fpr, tpr)
        #rounded = [round(x[0]) for x in predictions]
        #accuracy = accuracy_score(Y_cv_test,rounded)
        cv_outcomes_acc.append(accuracy)
        cv_outcomes_sep.append(sepcificity)
        cv_outcomes_sen.append(sensitivity)
        cv_outcomes_str.append(strength)
        cv_outcomes_mcc.append(mcc)
        cv_outcomes_f1.append(f1)
        cv_outcomes_auc.append(aucvalue)
        fileout.write("Cross_validation:\n")
        fileout.write("Accuracy: "+str(accuracy))
        fileout.write("\n")
        fileout.write("Sepcificity: "+str(sepcificity))
        fileout.write("\n")
        fileout.write("Sensitivity: "+str(sensitivity))
        fileout.write("\n")
        fileout.write("Strength: "+str(strength))
        fileout.write("\n")
        fileout.write("MCC: "+str(mcc))
        fileout.write("\n")
        fileout.write("F-score: "+str(f1))
        fileout.write("\n")
        fileout.write("AUC: "+str(aucvalue))
        fileout.write("\n")
        fileout.write("confusion matrix:")
        #fileout.write(str(confusion_matrix(Y_cv_test, rounded)))
        fileout.write(str(confusion_matrix(Y_cv_test, predictions)))
        fileout.write("\n")
        fileout.write("Classification report:")
        #fileout.write(str(classification_report(Y_cv_test, rounded)))
        fileout.write(str(classification_report(Y_cv_test, predictions)))
        fileout.write("\n")
    cv_acc = np.mean(cv_outcomes_acc)
    cv_acc_std = np.std(cv_outcomes_acc)
    cv_sep = np.mean(cv_outcomes_sep)
    cv_sen = np.mean(cv_outcomes_sen)
    cv_str = np.mean(cv_outcomes_str)
    cv_mcc = np.mean(cv_outcomes_mcc)
    cv_f1 = np.mean(cv_outcomes_f1)
    cv_auc = np.mean(cv_outcomes_auc)
    fileout.write("Accuracy_mean: "+str(cv_acc))
    fileout.write("\n")
    fileout.write("Accuracy_std: "+str(cv_acc_std))
    fileout.write("\n")
    fileout.write("Sepcificity_mean: "+str(cv_sep))
    fileout.write("\n")
    fileout.write("Sensitivity_mean: "+str(cv_sen))
    fileout.write("\n")
    fileout.write("Strength_mean: "+str(cv_str))
    fileout.write("\n")
    fileout.write("MCC_mean: "+str(cv_mcc))
    fileout.write("\n")
    fileout.write("Fscore_mean: "+str(cv_f1))
    fileout.write("\n")
    fileout.write("AUC_mean: "+str(cv_auc))
    fileout.write("\n")
    fileout.write("\n")
    fileout.write("\n")
    return cv_preds, cv_probs
        
def plotROC(Y_train, pred_train_prob, name, color):
    fpr, tpr, thresholds = roc_curve(Y_train, pred_train_prob)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, color=color, lw=1, label=name+' (area = %0.4f)' % roc_auc)
    plt.legend(loc='best')
    #plt.plot([0, 1], [0, 1],'r--', label = 'random guessing = 0.5')
    #plt.legend(loc = 'lower right')
    plot_name = 'ROC_'+name+'_AutoEncoder.png'
    plt.savefig(plot_name)

def machine_learning(X_pre,Y_train,random_seed,fileout):
    models = []
    #X_train=deep_autoencoder(X_pre)
    # raw data with no using autoencoder:
    X_train = X_pre
    para_lr = {'C' : [0.2,0.3,0.4,0.5,0.6], 
           'penalty' : ['l2'],
           'class_weight' : [{1:3.1471136}]
          }
    para_svm = {'C' : [128, 256, 512], 
            'kernel' : ['rbf'],
            'gamma' : [0.001,0.005,0.06324,0.09125,0.132633,0.3],
            'class_weight' : [{1:3.1471136}],
            'probability' : [True]
           }
    para_rf = {'n_estimators': [65,102,256,500], 
           'max_features': ['log2', 'sqrt'], 
           'criterion': ['entropy', 'gini'],
           'max_depth': [5,7,9,12], 
           'min_samples_split': [2,4,6,7],
           'class_weight' : [{1:3.1471136}]
          }
    
    para_mlp = {'alpha' : [0.01,0.1,1.0],
             'activation' : ['relu'],
             'max_iter' : [64,128],
             'hidden_layer_sizes' : [(256,128,32),(64,32,8),(128,64,16)],
             'solver' : ['lbfgs','sgd','adam'],
             'random_state' : [int(random_seed)]
            }
    
    best_para_lr = para_tuning(LogisticRegression(solver = 'lbfgs', max_iter=256), para_lr, X_train, Y_train)
    print("Parameter for logistic regression finished")
    print (best_para_lr)
    best_para_svm = para_tuning(SVC(), para_svm, X_train, Y_train)
    print("Parameter for svm finished")
    print (best_para_svm)
    best_para_rf = para_tuning(RandomForestClassifier(), para_rf, X_train, Y_train)
    print("Parameter for random forest finished")
    print (best_para_rf)
    best_para_mlp = para_tuning(MLPClassifier(), para_mlp, X_train, Y_train)
    print("Parameter for MLP finished")
    print (best_para_mlp)
    
    models.append(('LR', best_para_lr,'red'))
    models.append(('SVM', best_para_svm,'blue'))
    models.append(('RandomForest', best_para_rf,'green'))
    models.append(('ANN', best_para_mlp,'darkorange'))

    seed = int(random_seed)
    for name, model,color in models:
        fileout.write(name+": ")
        fileout.write(str(model))
        fileout.write("\n")
        cv_preds, cv_probs = cv(model, X_train, Y_train, seed, fileout)
        plotROC(np.asarray(cv_preds), np.asarray(cv_probs), name, color)
        '''model.fit(X_train, Y_train)
        predictions = model.predict(X_test)
        pred_train_prob = model.predict_proba(X_test)[:, 1]
        fileout.write("Prediction:\n")
        fileout.write("Accuracy: "+str(accuracy_score(Y_test, predictions)))
        fileout.write("\n")
        confusion = confusion_matrix(Y_test, predictions)
        TP = confusion[1, 1]
        TN = confusion[0, 0]
        FP = confusion[0, 1]
        FN = confusion[1, 0]
        sepcificity = TN / float( TN + FP)
        sensitivity = TP / float(FN + TP)
        strength = (sepcificity + sensitivity)/2
        mcc = matthews_corrcoef(Y_test, predictions)
        f1 = f1_score(Y_test, predictions)
        fpr, tpr, thresholds = roc_curve(Y_test, pred_train_prob)
        aucvalue = auc(fpr, tpr)
        fileout.write("Sepcificity: "+str(sepcificity))
        fileout.write("\n")
        fileout.write("Sensitivity: "+str(sensitivity))
        fileout.write("\n")
        fileout.write("Strength: "+str(strength))
        fileout.write("\n")
        fileout.write("MCC: "+str(mcc))
        fileout.write("\n")
        fileout.write("Fscore: "+str(f1))
        fileout.write("\n")
        fileout.write("AUC: "+str(aucvalue))
        fileout.write("\n")
        fileout.write("confusion matrix:")
        fileout.write(str(confusion_matrix(Y_test, predictions)))
        fileout.write("\n")
        fileout.write("Classification report:")
        fileout.write(str(classification_report(Y_test, predictions)))
        fileout.write("\n")
        plotROC(Y_test, pred_train_prob, name, color) 
        '''

def main():
    from sklearn.utils import shuffle 

    random_seed = 345
    fileout = open("models_parameter_tuning.txt", "w")
    ### log(RPKM+1) normalization expression data:
    exp_feature  = pd.read_csv("2_expression_normalized.csv", header = None)
    min_max_scaler = preprocessing.MinMaxScaler()
    exp_feature = min_max_scaler.fit_transform(exp_feature)
    exp_feature = pd.DataFrame(exp_feature)
    exp_label = pd.read_csv("2_labels.csv", header = None)
    exp_label.columns = ['label']

    exp = pd.concat([exp_feature, exp_label], axis=1)
    exp = shuffle(exp)
    
    # try embedding:
    exp = pd.read_csv("ASD_genes_embeded.csv")
    my_dict = {'ASD': 1, 'Disease': -1}
    labels = [my_dict[i] for i in exp['Gene.Type']]
    labels = pd.DataFrame(labels)
    labels.columns=['label']

    exp_f = exp.iloc[:,2:52] 
    min_max_scaler = preprocessing.MinMaxScaler()
    exp_f = min_max_scaler.fit_transform(exp_f)
    exp_f = pd.DataFrame(exp_f)

    exp = pd.concat([exp_f,labels],axis=1)
    exp = shuffle(exp,random_state=0)
    exp.index = range(2227)
    
    #
    Xdata_train = exp.drop('label', axis = 1)
    Xdata_train = Xdata_train.reset_index(drop = True)
    Xdata_train = np.array(Xdata_train)
    
    Ydata_train = exp['label'].reset_index(drop = True)
    Ydata_train = np.array(Ydata_train)
    


    machine_learning(Xdata_train,Ydata_train,random_seed,fileout)
    
if __name__ == '__main__':
        main()