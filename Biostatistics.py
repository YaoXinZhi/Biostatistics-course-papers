# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 07/06/2020 下午2:09
@Author: xinzhi yao
"""

import random
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.svm import SVC
from collections import defaultdict
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score

def read_fasta_file(classification_file: str, pos=False):
    seq_list = []
    with open(classification_file) as f:
        for line in f:
            l = line.strip()
            if not l.startswith('>') and l != '':
                seq_list.append(l.upper())

    if pos:
        labels = ['pos']*len(seq_list)
    else:
        labels = ['neg']*len(seq_list)
    print('sequence count: {0}.'.format(len(seq_list)))
    return seq_list, labels


def read_regression_file(regression_file: str):
    pos2y = {}
    with open(regression_file) as f:
        for line in f:
            l = line.strip().split('\t')
            pos2y[(l[1], (int(l[2]), int(l[3])))] = int(l[4])
    return pos2y


def cut_seq(seq, k):
    result = []
    for i in range(len(seq)-k+1):
        result.append(seq[i: i+k])
    return result


def kmer(seq: str, k: int):
    kmer_seq = [''.join(i) for i in itertools.product('ATCG', repeat=k)]
    value = ['0']*len(kmer_seq)
    kmer_count = dict(zip(kmer_seq, value))
    seq_cut = cut_seq(seq, k)
    for kmer in kmer_seq:
        num = seq_cut.count(kmer)
        kmer_count[kmer] = num
    total = sum(list(kmer_count.values()))
    kmer_feature = [freq/total for freq in list(kmer_count.values())]
    return kmer_feature


def make_dataset(seq_list):
    k = 6
    X = [kmer(seq, k) for seq in seq_list]
    X = np.array(X)
    return X


def classifation(model, parameter, X, Y):
    x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    clf = model()
    grid = GridSearchCV(clf, parameter, cv=10, n_jobs=4).fit(x_train, y_train)
    best = grid.best_estimator_
    predict = best.predict(x_test)
    acc = accuracy_score(y_test, predict)
    auc = roc_auc_score(y_test, predict)
    print('============={0}============='.format(model))
    print('best parameters: {0}'.format(grid.best_params_))
    print(classification_report(y_test, predict))
    print('acc: {0:.2f}'.format(acc))
    print('auc: {0:.2f}'.format(auc))


def add_label(input_file: str, out: str, label_num=113):
    wf = open(out, 'w')
    label_count = defaultdict(int)
    labels = [i for i in range(1, 11)]
    with open(input_file) as f:
        for line in f:
            l = line.strip().split(' ')
            label = random.sample(labels, 1)[0]
            label_count[label] += 1
            wf.write('{0},{1}\n'.format(','.join(l), label))
            if label_count[label] == label_num:
                labels.remove(label)
    wf.close()


def kmer_feature(seq_list: list, k: int, base: str, out: str):
    kmers = list(itertools.product(base, repeat=k))
    feature_matrix = []
    for i in range(len(kmers)):
        kmers[i] = ''.join(kmers[i])

    for seq in tqdm(seq_list):
        kmer_in_seq = [ ]
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            kmer_in_seq.append(kmer)
        count = []
        for i in kmers:
            count.append(kmer_in_seq.count(i))
        feature_vector = [str(c/len(kmer_in_seq)) for c in count]
        feature_matrix.append(feature_vector)

    feature_matrix = pd.DataFrame(feature_matrix)
    print(feature_matrix.shape)
    feature_matrix.to_csv(out, index=False, header=False, encoding='gbk', float_format='%.4f', sep="\t")
    print('save done: {0}'.format(out))

def feature_add_label(input_file: str, out: str, pos=True):
    if pos:
        label = 0
    else:
        label = 1
    wf = open(out, 'w')
    count = 1
    with open(input_file) as f:
        for line in f:
            l = line.strip().split()
            if count == 1 and pos:
                print(len(l))
                wf.write('{0},class\n'.format(','.join(['V'+str(i) for i in range(len(l))])))
                count += 1
            wf.write('{0},{1}\n'.format(','.join(l), label))
    wf.close()

def regression_add_Y(input_file: str, out: str, label_list: list):
    wf = open(out, 'w')
    count = 0
    with open(input_file) as f:
        for line in f:
            l = line.strip().split()
            if count == 0:
                wf.write('{0},Y\n'.format(','.join([ 'V' + str(i) for i in range(len(l))])))
            wf.write('{0},{1}\n'.format(','.join(l), label_list[count]))
            print(count)
            count += 1
    wf.close()

def read_label(input_file: str):
    label_list = []
    with open(input_file) as f:
        for line in f:
            l = line.strip().split()
            label_list.append(l[-1])
    print('label number: {0}.'.format(len(label_list)))
    return label_list

def metrics_confusion(TP: int, FN: int, FP: int, TN: int):
    accuracy = (TP+TN)/(TP+FN+TN+FP)
    precision = TP/(TP+FP)
    recall = TP/(TP+FN)
    F1_score = (precision*recall)/(precision+recall)
    print('Accuracy: {0:.5f}, Precision: {1:.5f}, Recall: {2:.5f}, F1-score: {3:.5f}'.\
          format(accuracy, precision, recall, F1_score))

def read_regression_result(input_file: str):
    value_list = []
    with open(input_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            value_list.append(float(l[-1]))
    print('value number: {0}.'.format(len(value_list)))
    return value_list

def metrics_regression(pred_file: str, ture_file: str):
    pred_list = read_regression_result(pred_file)
    true_list = read_regression_result(ture_file)

    MSE = mean_squared_error(true_list, pred_list)
    R2 = r2_score(true_list, pred_list)
    print('MSE: {0}, R2: {1}'.format(MSE, R2))



if __name__ == '__main__':

    classification_neg_file = './data/Classification/ABI5-neg.fasta'
    classification_pos_file = './data/Classification/ABI5-pos.fasta'
    classification_neg_feature_file = './data/Classification/ABI5-neg.feature.fasta'
    classification_pos_feature_file = './data/Classification/ABI5-pos.feature.fasta'
    classification_pos_feature_label_file = './data/Classification/ABI5-pos.feature.label.txt'
    classification_neg_feature_label_file = './data/Classification/ABI5-neg.feature.label.txt'

    regression_file = './data/Regression/p53-chipseq-1130sgRNA.txt'
    regression_file_label = './data/Regression/p53-chipseq-1130sgRNA.csv'
    regression_file_feature = './data/Regression/1130sgRNA_100bp_7mer.txt'
    regression_feature_file = './data/Regression/regression.feature.txt'
    # k number of kmer
    k = 6
    base = 'ATCG'
    example_fasta = './kx/Computational_CRISPR_Strategy-master/example/example.fasta'

    # 给数据打上标签 用于十折交叉验证
    add_label(regression_file, regression_file_label)
    # 通过kmer生成特征文件


    seq_list, _ = read_fasta_file(classification_pos_file)
    feature_matrix = kmer_feature(seq_list, k, base, classification_pos_feature_file)
    # seq_list, _ = read_fasta_file(classification_neg_file)
    # feature_matrix = kmer_feature(seq_list, k, base, classification_neg_feature_file)

    test_in = './data/Classification/test.txt'
    test_out = './data/Classification/test.label.txt'
    feature_add_label(test_in, test_out, True)

    feature_add_label(classification_neg_feature_file, classification_neg_feature_label_file, False)
    feature_add_label(classification_pos_feature_file, classification_pos_feature_label_file, True)

    # 给回归特征数据打上Y值标签
    label_list = read_label(regression_file)
    regression_add_Y(regression_file_feature, regression_feature_file, label_list)

    # pred_file = './result/lasso_pred/lasso.50.pred.txt'
    # ture_file = './result/lasso_pred/lasso.50.true.txt'
    # pred_file = './result/lasso_pred/lasso.100.pred.txt'
    # ture_file = './result/lasso_pred/lasso.100.true.txt'
    # pred_file = './result/lasso_pred/lasso.200.pred.txt'
    # ture_file = './result/lasso_pred/lasso.200.true.txt'

    pred_file = './result/ridge_pred/ridge.50.pred.txt'
    ture_file = './result/ridge_pred/ridge.50.true.txt'
    # pred_file = './result/ridge_pred/ridge.100.pred.txt'
    # ture_file = './result/ridge_pred/ridge.100.true.txt'
    # pred_file = './result/ridge_pred/ridge.200.pred.txt'
    # ture_file = './result/ridge_pred/ridge.200.true.txt'

    metrics_regression(pred_file, ture_file)

    # 读取数据
    # classification_neg_seq_list, label_neg = read_classification_file(classification_neg_file, pos=False)
    # classification_pos_seq_list, label_pos = read_classification_file(classification_pos_file, pos=True)
    # seq_list = classification_neg_seq_list + classification_pos_seq_list
    # label_list = label_neg + label_pos
    # X = make_dataset(seq_list)
    # Y = np.array(label_list)

    # 定义参数
    parameter_svm = {
        'kernel':('linear', 'rbf'),
        'C': [1, 2, 4, 10],
        'gamma': [0.125, 0.25, 0.5, 1, 2, 4]
    }

    parameter_tree = {
        'max_depth': range(1, 6)
    }

    parameter_rf = {
        'n_estimators': [3, 10, 15, 20],
        'criterion': ['gini', 'entropy'],
        'min_samples_leaf': [2, 4, 6],
    }

    parameter_lr = {
        'penalty': ['l1', 'l2'],
        'C': [0.01, 0.1, 1, 10]
    }

    parameter_nb = {
        'alpha': [1e-1, 1e-2]
    }

    # 模型训练与评价： 利用网格搜索方法进行寻参， 然后用测试集数据在最优模型进行测试，得到模型评价结果。
    model = [MultinomialNB, DecisionTreeClassifier, SVC, RandomForestClassifier, LogisticRegression]
    parameter = [parameter_nb, parameter_tree, parameter_svm, parameter_rf, parameter_lr]

    # for i, j in zip(model, parameter):
    #     classifation(i, j, X, Y)

    # 读取回归数据


