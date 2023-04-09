import os
import sys
import math
import random
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from pandas.core.frame import DataFrame

from matplotlib import pyplot as plt
from scipy.cluster import hierarchy
import scanpy as sc

from sklearn import svm
import joblib
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import AdaBoostClassifier, ExtraTreesClassifier, RandomForestClassifier
from sklearn.feature_selection import VarianceThreshold, f_classif, mutual_info_classif, SelectFromModel, SelectKBest
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix, auc, roc_curve, average_precision_score, f1_score, precision_recall_curve
from sklearn.calibration import CalibratedClassifierCV


def load_data(EXP_path, label_path):
    EXP = pd.read_csv(EXP_path, sep="\t", index_col=0)
    EXP = EXP.fillna(0)
    label = pd.read_csv(label_path, sep="\t", index_col=0)
    # print(EXP.head())
    # print(label.head())
    print("load data from " + EXP_path)
    print("sample_num:" + str(len(label)) + " " + "feature_num:" + str(len(EXP)))
    print("samlpe_classes:" + str(label["Phenotype"].unique()))  # ['PE' 'C']
    # 判断基因是否存在重复, 去除重复的基因, 不存在重复
    # print(len(GSE_EXP.index) == len(set(GSE_EXP.index)))  # TRUE
    # EXP = EXP[(EXP.T != 0).any()]  # 去除表达谱中全为零的行
    return EXP.T, label

# 基于阈值的梯度进行过滤
# def filter(GSE, GSE_number, threshold):
#     """ 基于F检验的有监督特征过滤
#     :return: 过滤后的选择的特征
#     """
#     f_values, f_pvalues = f_classif(GSE, GSE_number)
#     k = f_values.shape[0] - (f_pvalues > threshold).sum()
#     max = np.nanmax(f_pvalues)  # 忽略nan找到最大值
#     min = np.nanmin(f_pvalues)
#     if (threshold >= max) | (threshold <= min):
#         return DataFrame()
#     index = np.where(f_pvalues < threshold)[0].tolist()
#     return GSE[np.array(GSE.columns.tolist())[np.array(index)]]


# 基于基因数的阈值过滤
def filter(GSE, GSE_number, threshold):
    """ 基于F检验的有监督特征过滤
    :return: 过滤后的选择的特征
    """
    f_values, f_pvalues = f_classif(GSE, GSE_number)
    selector = SelectKBest(f_classif, k=4)
    selector.fit(GSE, GSE_number)

    dict_threshold = {}
    for i, j in zip(GSE.columns.values, f_pvalues):
        if str(j) == "nan":
            continue
        else:
            dict_threshold[i] = j
    dict_threshold = dict(sorted(dict_threshold.items(), key=lambda x: x[1], reverse=False))
    max = len(dict_threshold)
    min = 0
    if (threshold >= max) | (threshold <= min):
        return DataFrame(), 0
    return GSE[list(dict_threshold.keys())[:threshold]], list(dict_threshold.values())[threshold-1]


def vt_filter(data, threshold=0.0, number=None):
    """
    :param data: 表达谱数据框，行为样本，列为特征
    :param threshold: 最小方差阈值
    :return: 基于方差过滤后的数据框，及其全部基因的方差值
    """
    vt = VarianceThreshold(threshold=threshold)
    vt.fit(data)
    # 查看各个特征的方差
    dict_variance = {}
    for i, j in zip(data.columns.values, vt.variances_):
        dict_variance[i] = j
    # 获取保留了的特征的特征名
    ls = list()
    variances = list()
    flag = 0
    dict_variance = dict(sorted(dict_variance.items(), key=lambda x: x[1], reverse=True))
    # plt.hist(x=list(dict_variance.values()), bins=100, color='steelblue', edgecolor='black')
    # plt.xlabel('variance')
    # plt.ylabel("frequence")
    # plt.title('Variance distribution of genes')
    # plt.show()
    for i, j in dict_variance.items():
        if number != None:
            if flag < number:
                ls.append(i)
                variances.append(j)
                flag = flag + 1
        else:
            if j > threshold:
                ls.append(i)
                variances.append(j)
    # filter_data = pd.DataFrame(vt.fit_transform(data), columns=ls, index=data.index)
    filter_data = pd.DataFrame(data[ls], columns=ls, index=data.index)
    print("filter feature base on VarianceThreshold:" + str(threshold))
    print("inputshape:" + str(data.shape) + " " + "outputshape:" + str(filter_data.shape))

    return filter_data, variances


def model_filter(fitmodel, data):
    model = SelectFromModel(fitmodel, prefit=True)
    flags = model.get_support()
    select_feature = []
    for flag, feature in zip(flags, data.columns.tolist()):
        if flag == True:
            select_feature.append(feature)
    new_data = model.transform(data)
    new_data = DataFrame(new_data, index=data.index.tolist(), columns=select_feature)  # 数组变成数据框
    print("filter feature base on model:" + str(data.shape) + " to " + str(new_data.shape))
    return new_data


def cross_vaildation_svm(train_data, train_label, test_data, test_label):
    SVM_model = svm.SVC(random_state=1)
    # param_grid_linear = [{'kernel': ['linear'], 'C': np.arange(0.0, 0.002, 0.00002)}]
    param_grid_linear = [{'kernel': ['linear'], 'C': np.arange(0.0, 1.0, 0.02)}]

    grid = GridSearchCV(estimator=SVM_model, param_grid=param_grid_linear, scoring=None, n_jobs=8, cv=5)  # cv=None=3-fold cv

    grid.fit(train_data, train_label)

    best_params = list(grid.best_params_.items())
    print(best_params)
    towritedata = pd.DataFrame(grid.cv_results_)
    towritedata.to_csv('./modelselect/SVM_CVResults.csv')

    SVM_model1 = svm.SVC(C=1.0, kernel='rbf', random_state=1)
    SVM_model1.set_params(**grid.best_params_)
    SVM_model1.fit(train_data, train_label)
    cal = CalibratedClassifierCV(SVM_model1, cv=5)
    cal.fit(train_data, train_label)

    importance = SVM_model1.coef_
    joblib.dump(SVM_model1, './modelselect/SVM.sav')  # 保存模型q
    # 训练集上的分类指标
    y_predict_prob = cal.predict_proba(train_data)
    # Compute ROC curve and ROC area for each class
    fpr, tpr, thresholds = roc_curve(train_label, y_predict_prob[:, 1])
    roc_auc = auc(fpr, tpr)
    precision, recall, thresholds = precision_recall_curve(train_label, y_predict_prob[:, 1])
    f1 = 2 * (precision * recall) / (precision + recall)
    index = np.nanargmax(f1)
    train_index = [train_data.shape[1], f1[index], precision[index], recall[index], roc_auc, best_params]
    # 测试集上的分类指标
    y_predict_prob = cal.predict_proba(test_data)
    # Compute ROC curve and ROC area for each class
    fpr, tpr, thresholds = roc_curve(test_label, y_predict_prob[:, 1])
    roc_auc = auc(fpr, tpr)
    precision, recall, thresholds = precision_recall_curve(test_label, y_predict_prob[:, 1])
    f1 = 2*(precision*recall)/(precision+recall)
    index = np.nanargmax(f1)
    eva_index = [f1[index], precision[index], recall[index], roc_auc]
    # print("训练集：", SVM_model1.score(train_data, train_label))
    # print("测试集：", SVM_model1.score(test_data, test_label))
    return best_params, SVM_model1, importance, train_index, eva_index


def cross_vaildation_RF(train_data, train_label, test_data, test_label):
    # param_grid = {'n_estimators': range(1, 101, 10), "max_features": range(1, 11, 1)}
    param_grid = {'n_estimators': range(1, 101, 1), "max_features": range(1, 21, 1)}  # 0.94, 0.65
    # param_grid = {'n_estimators': range(1, 101, 1), "max_features": range(100, 201, 10)}  #
    clf0 = RandomForestClassifier(oob_score=False, random_state=1)
    grid_clf0 = GridSearchCV(clf0, param_grid, cv=5, n_jobs=10, return_train_score=True)
    grid_clf0.fit(train_data, train_label)

    towritedata = pd.DataFrame(grid_clf0.cv_results_)
    towritedata.to_csv('./modelselect/RF_CVResults1.csv')

    clf = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
                                 max_depth=2, max_features='auto', max_leaf_nodes=None,
                                 min_impurity_decrease=0.0, random_state=1,
                                 min_samples_leaf=1, min_samples_split=2,
                                 min_weight_fraction_leaf=0.0, n_estimators=20, n_jobs=1,
                                 oob_score=False, verbose=0, warm_start=False)
    best_param = list(grid_clf0.best_params_.items())
    print(best_param)
    clf.set_params(**grid_clf0.best_params_)
    clf.fit(train_data, train_label)
    importance = clf.feature_importances_
    print("训练集：", clf.score(train_data, train_label))
    print("测试集：", clf.score(test_data, test_label))
    print(classification_report(clf.predict(test_data), test_label))
    joblib.dump(clf, './modelselect/RF1.sav')
    n_classes = 2
    y_predict_prob = clf.predict_proba(test_data)

    # Compute ROC curve and ROC area for each class
    fpr, tpr, thresholds = roc_curve(test_label, y_predict_prob[:, 1])
    roc_auc = auc(fpr, tpr)
    print(roc_auc)

    return best_param[0][1], best_param[1][1], clf, importance, roc_auc


def evaluate(label, predict, visual=False):
    threshold = 0
    TPR_list = []
    FPR_list = []
    TY_list = []
    F1_score = []
    Recall_score = []
    precesion = []
    while threshold <= 1:  # 设置阈值
        y_predicted = []
        for i in range(0, len(label)):
            if predict[i] >= threshold:
                y_predicted.append(1)
            else:
                y_predicted.append(0)
        cm = confusion_matrix(label, y_predicted)
        TP = cm[1, 1]
        TN = cm[0, 0]
        FP = cm[0, 1]
        FN = cm[1, 0]
        FPR = round(FP / (FP + TN), 3)
        TPR = round(TP / (TP + FN), 3)
        TY = round(TN / (TN + FP), 3)
        prece = round(TP / (TP + FP), 3)
        TPR_list.append(TPR)
        FPR_list.append(FPR)
        TY_list.append(TY)
        precesion.append(prece)
        Recall_score.append(TPR)
        F1_score.append(2 * prece * TPR / (prece + TPR))
        threshold += 0.05
    precesion[-1] = 1
    F1_score[-1] = 2 * precesion[-1] * Recall_score[-1] / (precesion[-1] + Recall_score[-1])
    index = np.nanargmax(F1_score)
    info = str(round(F1_score[index], 3)) + '(' + str(FPR_list[index]) + ',' + str(TPR_list[index]) + ')'
    roc_auc = auc(FPR_list, TPR_list)

    plt.figure()
    lw = 2
    plt.figure(figsize=(7, 5))
    plt.plot(FPR_list, TPR_list, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([-0.05, 1.0])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic ')
    plt.legend(loc="lower right")
    plt.show()

    return F1_score[index], FPR_list[index], TPR_list[index], roc_auc


def main(labels, EXPs, thresholdes):
    # load data
    random.seed(1)
    num_label = [1 if label == "PE" else 0 for label in labels["Phenotype"]]
    # 统计不同阈值下过滤得到的基因数目
    SVM_eva_list, RF_eva_list, GNB_eva_list, use_threshold = [], [], [], []
    EXPs, variances = vt_filter(EXPs, threshold=0.0)   # feature select based on variance 无监督
    # for threshold in np.arange(0.0, 0.2, 0.01):
    for threshold in thresholdes:
    # for number in range(1000, 1001, 100):
        # 记录高可变基因数目对模型的影响
        # 有监督的挑选高可变基因, 利用阈值
        # filter_GSE = filter(T_GSE_EXP, GSE_number, threshold=threshold)
        # 有监督的挑选高可变基因, 利用特征数
        filter_exp, top_thresh = filter(EXPs, num_label, threshold=threshold)
        if filter_exp.empty:
            continue
        else:
            use_threshold.append(threshold)
            # filter_GSE, variances = vt_filter(T_GSE_EXP, number=number)
            # print("get high variance gene top " + str(number))

            # region 可视化
            # obs = DataFrame(index=filter_GSE.index.tolist())
            # var = DataFrame(index=filter_GSE.columns.tolist())
            # adata = sc.AnnData(filter_GSE, obs=obs, var=var)  # 经过方差过滤掉的基因
            # sc.pp.scale(adata, max_value=10)  # 数据缩放，为了放到pca中进行可视化
            # adata.obs["catalog"] = GSE_label  # 添加类别标签
            # adata.obs["spec_catalog"] = specific_GSE_label  # 添加更细的类别标签
            # sc.tl.pca(adata, svd_solver='arpack')
            # sc.pl.pca(adata, color='catalog')
            # # sc.pl.pca_variance_ratio(adata, log=True)
            # # adata.write(results_file)
            # # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
            # # sc.tl.umap(adata)
            # # sc.pl.umap(adata, color=['catalog'])
            # # # sc.pl.umap(adata, color=['ABAT', 'ABCA12', 'ABCA4'])
            # # # sc.pl.umap(adata, color=['ABAT', 'ABCA12', 'ABCA4'], use_raw=False)
            # # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            # # adata.raw = adata  # 保存一下原始数据
            # 层次聚类
            # GSE_link = hierarchy.linkage(T_GSE_EXP, method='ward', metric='euclidean')
            # GSE_dend = hierarchy.dendrogram(GSE_link, labels=GSE_label["Phenotype"].tolist())
            # plt.show()
            # # 确定阈值得到聚类的结果
            # label = hierarchy.cut_tree(GSE_link, height=220)
            # label = label.reshape(label.size, )
            # PMID_link = hierarchy.linkage(T_PMID_EXP, method='ward', metric='euclidean')
            # PMID_dend = hierarchy.dendrogram(PMID_link, labels=PMID_label["Phenotype"].tolist())
            # plt.show()
            # endregion

            # split train and test
            train_data, test_data, train_label, test_label = train_test_split(
                filter_exp, num_label, random_state=1, train_size=0.7, test_size=0.3)

            # svm train
            best_params, best_model_svm, importance, train_index, eva_index = cross_vaildation_svm(train_data, train_label, test_data, test_label)
            SVM_eva_list.append(train_index + eva_index)
            # 从网格搜索得到的最佳模型中选出的高可变基因
            # new_GSE = model_filter(best_model_svm, filter_GSE)
            # GSE_link_new = hierarchy.linkage(new_GSE, method='ward', metric='euclidean')
            # GSE_dend = hierarchy.dendrogram(GSE_link_new, labels=GSE_label)
            # plt.show()

            # RF train
            # best_para1, best_para2, best_model_RF, importance, roc = cross_vaildation_RF(train_data, train_label, test_data, test_label)

            # region 决策树
            # clf = ExtraTreesClassifier()
            # ctree = clf.fit(train_data, train_label)
            # print("训练集：", clf.score(train_data, train_label))
            # print("测试集：", clf.score(test_data, test_label))
            # print(ctree.feature_importances_)
            #
            # model = SelectFromModel(ctree, prefit=True)
            # flags = model.get_support()
            # select_feature = []
            # for flag, feature in zip(flags, filter_GSE.columns.tolist()):
            #     if flag == True:
            #         select_feature.append(feature)
            # new_GSE = model.transform(filter_GSE)
            # new_GSE = DataFrame(new_GSE, index=filter_GSE.index.tolist(), columns=select_feature)  # 数组变成数据框
            # print(new_GSE.shape)
            #
            # obs = DataFrame(index=new_GSE.index.tolist())
            # var = DataFrame(index=new_GSE.columns.tolist())
            # adata = sc.AnnData(new_GSE, obs=obs, var=var)
            # adata.obs["catalog"] = GSE_label
            # sc.pp.scale(adata, max_value=10)
            # sc.tl.pca(adata, svd_solver='arpack')
            # sc.pl.pca(adata, color='catalog')
            # sc.pl.pca_variance_ratio(adata, log=True)
            # adata.write(results_file)
            #
            # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
            # sc.tl.umap(adata)
            # sc.pl.umap(adata, color=['catalog'])
            #
            # GSE_link_new = hierarchy.linkage(new_GSE, method='ward', metric='euclidean')
            # GSE_dend = hierarchy.dendrogram(GSE_link_new, labels=GSE_label)
            # plt.show()
            # endregion

            # GaussianNB train
            cls = GaussianNB()
            cls.fit(train_data, train_label)
            # 训练集上的分类指标
            y_predict_prob = cls.predict_proba(train_data)
            # Compute ROC curve and ROC area for each class
            fpr, tpr, thresholds = roc_curve(train_label, y_predict_prob[:, 1])
            roc_auc = auc(fpr, tpr)
            precision, recall, thresholds = precision_recall_curve(train_label, y_predict_prob[:, 1])
            f1 = 2 * (precision * recall) / (precision + recall)
            index = np.nanargmax(f1)
            train_index = [train_data.shape[1], f1[index], precision[index], recall[index], roc_auc]
            # 测试集上的分类指标
            y_predict_prob = cls.predict_proba(test_data)
            # Compute ROC curve and ROC area for each class
            fpr, tpr, thresholds = roc_curve(test_label, y_predict_prob[:, 1])
            roc_auc = auc(fpr, tpr)
            precision, recall, thresholds = precision_recall_curve(test_label, y_predict_prob[:, 1])
            f1 = 2 * (precision * recall) / (precision + recall)
            index = np.nanargmax(f1)
            eva_index = [f1[index], precision[index], recall[index], roc_auc]
            GNB_eva_list.append(train_index + eva_index)
            ## 返回经过模型挑选出的基因的权重
            # specific_GSE_path = "./data/bulk/GSE75010/ga_META_group_GSE.txt"
            # specific_GSE_label = pd.read_csv(specific_GSE_path, sep="\t", index_col=0)  # PE, C, LO
            # y_pred_train = cls.predict(train_data)
            # y_pred_test = cls.predict(test_data)
            # predict = y_pred_train
            # flagz, flagw, flagLO, flagPE = 0, 0, 0, 0
            # for label, pre, name in zip(train_label, predict, train_data.index.tolist()):
            #     s_label = specific_GSE_label.loc[name, "Phenotype"]
            #     if s_label == "LO":
            #         flagLO = flagLO + 1
            #     if s_label == "PE":
            #         flagPE = flagPE + 1
            #     if (pre != label) & (s_label == "LO"):
            #         flagw = flagw + 1
            #     if (pre != label) & (s_label == "PE"):
            #         flagz = flagz + 1
            # print(flagz, flagw, flagPE, flagLO)
            #
            # predict = y_pred_test
            # flagz, flagw, flagLO, flagPE = 0, 0, 0, 0
            # for label, pre, name in zip(test_label, predict, test_data.index.tolist()):
            #     s_label = specific_GSE_label.loc[name, "Phenotype"]
            #     if s_label == "LO":
            #         flagLO = flagLO + 1
            #     if s_label == "PE":
            #         flagPE = flagPE + 1
            #     if (pre != label) & (s_label == "LO"):
            #         flagw = flagw + 1
            #     if (pre != label) & (s_label == "PE"):
            #         flagz = flagz + 1
            # print(flagz, flagw, flagPE, flagLO)
    SVM_eva = DataFrame(SVM_eva_list, index=use_threshold, columns=["feature", "train_f1", "train_precision",
                                                                 "train_recall", "train_roc_auc", "best_para",
                                                                 "test_f1", "test_precision", "test_recall",
                                                                 "test_roc_auc"])
    GNB_eva = DataFrame(GNB_eva_list, index=use_threshold, columns=["feature", "train_f1", "train_precision",
                                                                 "train_recall", "train_roc_auc", "test_f1",
                                                                 "test_precision", "test_recall", "test_roc_auc"])
    return SVM_eva, GNB_eva


def is_image_file(filename):
    return any(filename.endswith(extension) for extension in ['.png', '.PNG', '.JPG', '.JPEG'])


if __name__ == '__main__':
    # input_path = sys.argv[1]
    # result_path = sys.argv[2]
    input_path = "PB_sig_150"
    result_path = "modelselect_PB_sig_150"

    GSE_EXP_path, GSE_label_path = "./data/bulk/GSE75010/EXP.txt", "./data/bulk/GSE75010/PAT_group.txt"
    PMID_EXP_path, PMID_label_path = "./data/bulk/PMID_33754042/EXP.txt", "./data/bulk/PMID_33754042/PAT_group.txt"
    GSE_EXP, GSE_label = load_data(GSE_EXP_path, GSE_label_path)
    PMID_EXP, PMID_label = load_data(PMID_EXP_path, PMID_label_path)
    # specific_GSE_path = "./data/bulk/GSE75010/ga_META_group_GSE.txt"
    # specific_GSE_label = pd.read_csv(specific_GSE_path, sep="\t", index_col=0)["Phenotype"].tolist()  # PE, C, LO
    print("GSE with PMID the same gene number:" + str(len(set(GSE_EXP.columns) & set(PMID_EXP.columns))))  # 14072
    print("="*50)
    print("train GSE bulk data without separate Z from W ")
    # thresholdes = range(100, 4001, 100)
    # SVM_eva, GNB_eva = main(GSE_label, GSE_EXP, thresholdes)
    # SVM_eva.to_csv("./" + result_path + "/SVM_eva.csv", sep="\t")
    # GNB_eva.to_csv("./" + result_path + "/GNB_eva.csv", sep="\t")

    print("=" * 50)
    print("train bulk data with separate Z from W ")
    zw_label = pd.read_csv("./data/bulk/GSE75010/zw_group_GSE.txt", sep="\t", index_col=0)
    z_label = zw_label[(zw_label["label"] == "Z_PE") | (zw_label["label"] == "Z_C")].index.tolist()
    w_label = zw_label[(zw_label["label"] == "W_PE") | (zw_label["label"] == "W_C")].index.tolist()
    Z_GSE, Z_label = GSE_EXP.loc[z_label], GSE_label.loc[z_label, "Phenotype"].tolist()
    W_GSE, W_label = GSE_EXP.loc[w_label], GSE_label.loc[GSE_EXP.loc[w_label].index.tolist()]["Phenotype"].tolist()

    # thresholdes = range(100, 4001, 100)
    # # thresholdes = np.append(np.arange(0.0, 0.1, 0.005), np.arange(0.1, 1, 0.05))
    # Z_SVM_eva, Z_GNB_eva = main(Z_label, Z_GSE, thresholdes)
    # Z_SVM_eva.to_csv("./" + result_path + "/Z_SVM_eva.csv", sep="\t")
    # Z_GNB_eva.to_csv("./" + result_path + "/Z_GNB_eva.csv", sep="\t")

    # W_SVM_eva, W_GNB_eva = main(W_label, W_GSE, thresholdes)
    # W_SVM_eva.to_csv("./" + result_path + "/W_SVM_eva.csv", sep="\t")
    # W_GNB_eva.to_csv("./" + result_path + "/W_GNB_eva.csv", sep="\t")

    ########################################################换用样本deconvolution表达的数据的结果
    single_dir = "./data/single/" + input_path
    cell_names = [x.split("_")[3] for x in os.listdir(single_dir) if is_image_file(x)]
    for cell_name in cell_names:
        path = "./data/single/" + input_path + "/CIBERSORTxHiRes_NA_" + cell_name + "_Window20.txt"
        # EXP = pd.read_csv(SCT_path, sep="\t", index_col=0)
        # print(len(np.where(np.isnan(EXP))[0]))
        # print(len(set(np.where(np.isnan(EXP))[0])))  # 1122个gene全为nan
        # EXP = EXP.fillna(0)
        print("=" * 50)
        print("single cell data result" + path)
        GSE_EXP, GSE_label = load_data(path, GSE_label_path)
        print(len(set(np.where(np.isnan(GSE_EXP))[1])))  # 打印数据中为nan的全部列
        obs = DataFrame(index=GSE_EXP.index.tolist())
        var = DataFrame(index=GSE_EXP.columns.tolist())
        adata = sc.AnnData(GSE_EXP, obs=obs, var=var)
        # adata.obs["catalog"] = GSE_label
        # adata.obs["spec_catalog"] = specific_GSE_label
        # sc.pp.normalize_total(adata, target_sum=1e4)  # 去文库的,bulk无需使用
        sc.pp.log1p(adata)
        # 挑选差异基因需要在scale之前做，去除掉通过方差过滤掉的基因
        sc_EXP = DataFrame(adata.X, index=GSE_EXP.index.tolist(),
                           columns=GSE_EXP.columns.tolist())  # 经过log的原始数据
        # NA_list = np.where(np.isnan(sc_EXP))[1]  # 找出为NA的列的索引
        # print(len(NA_list) / len(sc_EXP))
        # NA_list = NA_list[:int(len(NA_list) / len(sc_EXP))]
        # use_list = list(set(range(0, len(T_GSE_EXP.columns))).difference(set(NA_list)))
        sc_EXP = sc_EXP.dropna(axis=1, how='all')  # 去除全部为NA的列

        # main(GSE_EXP, GSE_label, sc_EXP)
        thresholdes = range(100, 4001, 100)
        SVM_eva, GNB_eva = main(GSE_label, sc_EXP, thresholdes)
        SVM_eva.to_csv("./" + result_path + "/" + cell_name + "_SVM_eva.csv", sep="\t")
        GNB_eva.to_csv("./" + result_path + "/" + cell_name + "_GNB_eva.csv", sep="\t")

        ##############################################################单细胞分早发晚发进行测试
        # Z_GSE, Z_label, T_Z_GSE = sc_EXP.loc[z_label], GSE_label.loc[sc_EXP.loc[z_label].index.tolist()][
        #     "Phenotype"].tolist(), sc_EXP.loc[z_label].T
        # W_GSE, W_label, T_W_GSE = sc_EXP.loc[w_label], GSE_label.loc[sc_EXP.loc[w_label].index.tolist()][
        #     "Phenotype"].tolist(), sc_EXP.loc[w_label].T

        # # thresholdes = np.append(np.arange(0.0, 0.1, 0.005), np.arange(0.1, 1, 0.05))
        # thresholdes = range(100, 4001, 100)
        # Z_SVM_eva, Z_GNB_eva = main(Z_label, Z_GSE, thresholdes)
        # Z_SVM_eva.to_csv("./" + result_path + "/Z_" + cell_name + "_SVM_eva.csv", sep="\t")
        # Z_GNB_eva.to_csv("./" + result_path + "/Z_" + cell_name + "_GNB_eva.csv", sep="\t")
        # # plt.plot(use_threshold, record, 'ro-', color='#4169E1', alpha=0.8, label='num_genes')
        # # plt.xlabel('EVT threhold')
        # # plt.ylabel("number of genes")
        # # plt.title('The number of genes that pass under the variance threshold')
        # # plt.show()

        # W_SVM_eva, W_GNB_eva = main(W_label, W_GSE, thresholdes)
        # W_SVM_eva.to_csv("./" + result_path + "/W_" + cell_name + "_SVM_eva.csv", sep="\t")
        # W_GNB_eva.to_csv("./" + result_path + "/W_" + cell_name + "_GNB_eva.csv", sep="\t")
        # # plt.plot(use_threshold, record, 'ro-', color='#4169E1', alpha=0.8, label='num_genes')
        # # plt.xlabel('EVT threhold')
        # # plt.ylabel("number of genes")
        # # plt.title('The number of genes that pass under the variance threshold')
        # # plt.show()












