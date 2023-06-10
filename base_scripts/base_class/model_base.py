#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-14

import os
import pandas as pd
import numpy as np
import shap
from base_scripts.base_class.global_base import model_control
from sklearn.feature_selection import VarianceThreshold, SelectFromModel, SelectKBest, chi2, f_classif, \
    mutual_info_classif
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler, MinMaxScaler, MaxAbsScaler, RobustScaler, Normalizer
from sklearn.metrics import f1_score, accuracy_score, roc_auc_score, confusion_matrix, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from hyperopt import hp, fmin, tpe, Trials, STATUS_OK
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator


class ml_model(model_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 从表单转换到对象的映射
        self.form2class = {
            # 缩放对象
            'StandardScaler': StandardScaler,
            'MinMaxScaler': MinMaxScaler,
            'MaxAbsScaler': MaxAbsScaler,
            'Normalizer': Normalizer,
            'RobustScaler': RobustScaler,
            # 特征选择器
            'SelectFromModel': SelectFromModel,
            'Chi-squared': chi2,
            'ANOVA': f_classif,
            'Mutual_info': mutual_info_classif,
            # 超参数寻优数据生成对象
            'uniform': hp.uniform,
            'uniformint': hp.uniformint,
            'quniform': hp.quniform,
            'loguniform': hp.loguniform,
            'qloguniform': hp.qloguniform,
            'randint': hp.randint,
            'choice': hp.choice,
        }

    def get_df(self, csv_file=None):
        # 若未指定CSV，则读取descriptor
        if not csv_file:
            csv_file = self.descriptor_csv
        # 读取描述符数据
        df = pd.read_csv(csv_file, encoding='utf-8')
        # 返回数据
        return df

    def get_x_y(self, df):
        # 根据df是否含有class列判断是否分成x,y
        if 'class' in df.columns:
            #  划分XY
            x = df.iloc[:, 1:-1]
            y = df.iloc[:, -1]
        else:
            x = df.iloc[:, 1:]
            y = None
        return x, y, np.array(x.columns)

    def train_test_split(self, *, x, y, train_size=0.8):
        # 判断是否分割数据集
        if train_size == 1:  # 若训练集数据为1，直接返回x,y
            return x, y
        else:  # 若训练集数比例不为1， 进行划分
            # 划分数据集
            train_x, test_x, train_y, test_y = train_test_split(x, y, train_size=train_size, shuffle=True,
                                                                stratify=y)  # 随机按照4：1划分训练集和测试集
            # 返回
            return train_x, test_x, train_y, test_y

    def data_scaler(self, scaler, *, train_x, test_x):
        # 实例化对象
        scaler = self.form2class[scaler]().fit(train_x)
        # 进行缩放
        train_x = scaler.transform(train_x)
        test_x = scaler.transform(test_x)
        return train_x, test_x

    def remove_same_feature(self, *, train_x, test_x, feature_names):
        # 设置阈值
        threshold = VarianceThreshold().fit(train_x)
        # 低方差过滤
        train_x = threshold.transform(train_x)
        test_x = threshold.transform(test_x)
        # 获取特征名
        feature_names = feature_names[threshold.get_support()]
        # 返回
        return train_x, test_x, feature_names

    def feature_selection(self, selector, *, train_x, train_y, test_x, para, feature_names):
        # 若选择不为None， 则继续执行
        if selector == 'SelectFromModel':  # 基于树的特征选择
            # 构建树
            clf = RandomForestClassifier(n_estimators=int(para), n_jobs=28)
            clf.fit(train_x, train_y)
            # 选择
            selector = self.form2class[selector](clf, prefit=True)
            train_x = selector.transform(train_x)
            test_x = selector.transform(test_x)
            # 获取特征名
            feature_names = feature_names[selector.get_support()]
        elif selector == 'Retain_all':  # 若为Retain_all，则不进行特征选择
            pass
        else:
            metric = self.form2class[selector]  # 统计检验指标
            selector = SelectKBest(metric, k=min(int(para), len(train_x))).fit(train_x, train_y)  # 构建选择器
            # 转换
            train_x = selector.transform(train_x)
            test_x = selector.transform(test_x)
            # 获取特征名
            feature_names = feature_names[selector.get_support()]
        return train_x, test_x, feature_names

    def cal_feature_importance(self, feature_names, x, y):
        # 实例化随机森林分类器
        clf = RandomForestClassifier(n_jobs=-1)
        # 训练数据
        clf.fit(x, y)
        # 输出重要性
        importance = clf.feature_importances_
        # 合并特征和重要性
        feature_importance = zip(feature_names, importance)
        # 排序
        feature_importance = sorted(feature_importance, key=lambda x: x[1], reverse=True)
        # 写入CSV
        pd.DataFrame(feature_importance).T.to_csv(self.feature_importance_csv, index=False, header=False, mode='a')

    # 判断数值是否在区间中, 若不在区间中，则返回默认值，在区间中则转换成合适的数字类型，返回
    def isin_range(self, value, default_value, para_range):
        # 判断是否为数字
        try:
            float(value)
        except:  # 非数字
            # 判断是否在区间内 max_features
            if value in para_range[2:]:
                return value
            else:
                return default_value
        else:  # 为数字
            # 区间取到0， 数值为浮点数
            if para_range[0] == 0:
                value = float(value)
                if 0 <= value <= para_range[1]:
                    return value
                else:
                    return default_value
            # 区间取不到0，且数值为浮点数
            elif para_range[0] == 0.1:
                value = float(value)
                if 0 < value <= para_range[1]:
                    return value
                else:
                    return default_value
            # 区间>=1，为整数
            else:
                value = int(value)
                if para_range[0] <= value <= para_range[1]:
                    return value
                else:
                    return default_value

    def check_range(self, value_list, model_default_value, para_range):
        # 区间数
        first_num = value_list[0]  # 第一个数
        last_num = value_list[-1]  # 第二个数
        first_num = self.isin_range(first_num, model_default_value, para_range)  # 第一个数
        last_num = self.isin_range(last_num, model_default_value, para_range)  # 第二个数
        minimum, maximum = min(first_num, last_num), max(first_num, last_num)  # 返回范围
        # 若参数值有问题，返回的都是默认值，则选用默认参数范围
        if minimum == maximum:
            minimum, maximum = para_range[0], para_range[1]
        return minimum, maximum

    def get_space(self, *, job_form, hyper_parameters, hyper_parameters_range):  # 获得超参数空间并验证超参数范围
        # 初始化space
        space = {}
        # 循环判断超参数是否需要寻优以及寻优方法
        for hyper, default_value in hyper_parameters.items():  # 获取超参数类型和超参数默认值
            value = job_form[hyper]  # 超参数值
            key = hyper  # 键 为超参数类型
            para_range = hyper_parameters_range.get(key)  # 参数合理范围
            if '-' in value:  # 若 - 存在，则进行寻优
                value_list = value.split('-')
                method = self.form2class[job_form[hyper + '_method']]  # 获取生成分布的方法 如 hp.uniform
                # 获取并检查区间范围
                minimum, maximum = self.check_range(value_list=value_list,
                                                    model_default_value=default_value,
                                                    para_range=para_range)
                # 根据生成分布方法不同，输入数据格式不同
                if method == hp.uniform or method == hp.uniformint:  # 若为均匀(整数)分布
                    space[key] = method(key, minimum, maximum)
                elif method == hp.loguniform:  # 若为指数均匀分布 最小值不能为0
                    # 判断是否为0
                    if minimum == 0:
                        minimum = 1*pow(10, -10)  # 浮点数最小值
                    minimum = np.log(minimum)
                    maximum = np.log(maximum)
                    space[key] = method(key, minimum, maximum)
                elif method == hp.randint:  # 随机整数 最后一位数取不到，所以+1 [low, high)
                    space[key] = method(key, minimum, maximum + 1)
                elif method == hp.choice:
                    space[key] = method(key, [int(i) if i.isdigit() else i for i in value_list])
            # 不进行寻优, 转换参数格式
            else:
                value = self.isin_range(value, default_value, para_range=para_range)  # 检查并合理化范围
                space[key] = value  # 赋值
        return space

    def hyperopt_tuning(self, *, train_x, train_y, model, parameter_space, max_evals):
        # 超参数寻优模型
        def tune_model(hyper_parameter):  # 待寻优函数
            try:  # xgb,rf有并行，若为svm，则报错
                clf = model(**hyper_parameter, random_state=42, n_jobs=-1)
            except:  # 为svm
                clf = model(**hyper_parameter, random_state=42)
            # 打印参数范围， 调试时使用
            # import hyperopt.pyll.stochastic
            # print(hyperopt.pyll.stochastic.sample(hyper_parameter))
            # 十折交叉验证  {'fit_time':[], 'score_time':[], 'test_accuracy':[], 'test_f1':[], 'train_accuracy':[],  'train_f1':[]}
            scores = cross_validate(clf, train_x, train_y, scoring=('f1', 'accuracy'),
                                    cv=StratifiedKFold(n_splits=10, shuffle=True, random_state=42), n_jobs=10,
                                    return_train_score=True)
            # 写数据到CSV中 训练集精度 验证集精度 验证集F1score
            validation_f1 = scores['test_f1'].mean()
            pd.DataFrame([scores['train_accuracy'].mean(), scores['test_accuracy'].mean(), validation_f1]).T.to_csv(
                self.bias_var_csv, index=False,
                header=False, mode='a')
            # 返回均值
            return -validation_f1

        # 实例化记录对象
        trials = Trials()
        # 寻优
        best_parameters = fmin(tune_model, parameter_space, algo=tpe.suggest, max_evals=max_evals,
                               rstate=np.random.RandomState(42), trials=trials)  # 寻找model()函数的最小值
        return best_parameters, trials

    # 方差偏差分析 作图
    def bias_var_analysis(self):
        # 读取数据
        df = pd.read_csv(self.bias_var_csv, encoding='utf-8', header=None)
        # 存成列表
        var = df.iloc[:, 0].tolist()  # 方差  在训练集上的精度
        bias = df.iloc[:, 1].tolist()  # 偏差  在测试集上的精度
        # 作图
        fig, ax = plt.subplots(figsize=(6, 6))  # 创建画布
        plt.plot(var, '-', label='Training_set')  # 方差线
        plt.plot(bias, '-', label='Validation_set')  # 偏差线
        ax.set(xlabel='Tuning Round', ylabel='Accuracy',  # 设置坐标轴标题
               title='Validation Curve')
        plt.legend()  # 设置标签
        # 控制坐标刻度为1
        # ax.xaxis.set_major_locator(MultipleLocator(1))
        # 自动控制间距
        plt.tight_layout()
        plt.savefig(self.bias_var_png, format='png', dpi=600, transparent=True)  # 储存图片
        # 关闭图片
        plt.close()

    # 做学习曲线
    def plot_learning_curve(self):
        # 读取数据
        df = pd.read_csv(self.bias_var_csv, encoding='utf-8', header=None)
        # 存成列表
        f1 = df.iloc[:, 2].tolist()  # f1score
        # 作图
        fig, ax = plt.subplots(figsize=(6, 6))  # 创建画布
        plt.plot(f1, '-')  # f1score
        ax.set(xlabel='Tuning Round', ylabel='F1_score',  # 设置坐标轴标题
               title='Learning Curve')
        # 控制坐标刻度为1
        # ax.xaxis.set_major_locator(MultipleLocator(1))
        # 自动控制间距
        plt.tight_layout()
        plt.savefig(self.learning_curve_png, format='png', dpi=600, transparent=True)  # 储存图片
        # 关闭图片
        plt.close()

    # 做超参数寻优图
    def plot_hyperParameter_scatter(self, trials):
        # 获取被寻优的超参数
        hypered_parameter = trials.trials[0]['misc']['vals']
        # 获取做的图的数量
        fig_num = len(hypered_parameter)
        # 获取彩虹色
        cmap = plt.cm.jet
        # 创建子图
        f, axes = plt.subplots(nrows=2, ncols=fig_num, figsize=(fig_num * 5, 10))
        # 循环创建图
        for i, val in enumerate(hypered_parameter):
            # 获取数据
            xs = np.array([t['misc']['vals'][val] for t in trials.trials]).ravel()  # 超参数取值
            ys = [-t['result']['loss'] for t in trials.trials]  # loss值
            zs = [t['tid'] for t in trials.trials]  # 时间步
            # 若调多个参数
            try:
                axes[0][i].scatter(zs, xs, s=20, linewidth=0.01, alpha=0.5,
                                   c=cmap(float(i) / len(hypered_parameter)))  # 画散点
                axes[1][i].scatter(xs, ys, s=20, linewidth=0.01, alpha=0.5,
                                   c=cmap(float(i) / len(hypered_parameter)))  # 画散点
                axes[0][i].set(xlabel='Tuning Round', ylabel='Value', title=val)  # 设置标题 参数值vs时间步
                axes[1][i].set(xlabel='Hyper_paramter_value', ylabel='-Loss', title=val)  # 设置标题  模型误差vs参数值
            # 调单个参数
            except:
                axes[0].scatter(zs, xs, s=20, linewidth=0.01, alpha=0.5,
                                c=cmap(float(i) / len(hypered_parameter)))  # 画散点
                axes[1].scatter(xs, ys, s=20, linewidth=0.01, alpha=0.5,
                                c=cmap(float(i) / len(hypered_parameter)))  # 画散点
                axes[0].set(xlabel='Tuning Round', ylabel='Value', title=val)  # 设置标题 参数值vs时间步
                axes[1].set(xlabel='Hyper_paramter_value', ylabel='-Loss', title=val)  # 设置标题  模型误差vs参数值
        # 自动控制间距
        plt.tight_layout()
        # 输出图片
        f.savefig(self.hyperParameter_scatter_png, format='png', dpi=600, transparent=True)
        # 关闭图片
        plt.close()

    def best2parameter(self, *, job_form, hyper_parameters, best_parameters):
        # 初始化最后参数字典
        fina_parameters = {}
        # 循环选取模型超参数
        for key in hyper_parameters:
            # 判断寻优参数中是否有该参数
            if key in best_parameters:  # 若有，则进一步判断其方法是否为choice
                if job_form[key + '_method'] == 'choice':  # 若是choice， 则生成对应列表，并赋值
                    value_list = [int(i) if i.isdigit() else i for i in job_form[key].split('-')]  # 可选择列表
                    fina_parameters[key] = value_list[best_parameters[key]]  # 根据最优参数选择
                else:  # 若不是choice，则赋值
                    fina_parameters[key] = best_parameters[key]
                pass
            else:  # 若无， 直接组成字典
                fina_parameters[key] = job_form[key]
        # 转换数据格式
        # fina_parameters = {i: float(j) if j.isdigit() else j for i, j in fina_parameters.items()}
        # 返回最终参数列表
        return fina_parameters

    # 计算富集因子
    def cal_ef(self, test_y, y_pred, *, top=0.01):
        # 合并df
        df = pd.DataFrame(test_y)  # 把真实值转成Dataframe
        df['pred'] = y_pred
        # 根据打分排序，并且给缺失值填充0
        df.sort_values(by='pred', inplace=True, ascending=False)
        # 去除同分异构分子，只保留第一项
        # data.drop_duplicates(subset='title', keep='first', inplace=True)
        # 总分子数
        N_total = len(df)
        # 总活性分子数
        N_active = len(df[df.iloc[:, 0] == 1])
        # 前b%总分子数
        topb_total = int(N_total * top + 0.5)
        # 前b%总分子数据
        topb_data = df.iloc[:topb_total, :]
        # 前b%总活性分子数
        topb_active = len(topb_data[topb_data.iloc[:, -1] >= 0.5])
        # 富集因子
        try:  # 防止因为数量太少导致乘上1%后为0
            ef_b = (topb_active / topb_total) / (N_active / N_total)
        except:
            ef_b = '-'
        return ef_b

    def model_metric(self, *, y_true, y_pred, y_pred_prob):
        # 把y_true 转成列表
        y_true = y_true.tolist()
        # 获取标签为1的概率 初始列表为 [ [0.6, 0.4], [0.4, 0.6], ...]
        y_pred_prob = [proba_for_zero_one[1] for proba_for_zero_one in y_pred_prob]
        # f1score
        f1 = f1_score(y_true=y_true, y_pred=y_pred)
        # 精度
        acc = accuracy_score(y_true=y_true, y_pred=y_pred)
        # roc_auc
        roc = roc_auc_score(y_true=y_true, y_score=y_pred_prob)
        # 混淆矩阵
        tn, fp, fn, tp = confusion_matrix(y_true=y_true, y_pred=y_pred).ravel()
        # 富集因子 1%
        ef = self.cal_ef(test_y=y_true, y_pred=y_pred_prob, top=0.01)
        # 标题栏
        columns = ['tn', 'fp', 'fn', 'tp', 'acc', 'roc', 'f1', 'ef1%', 'y_pred', 'y_pred_proba', 'y_true']
        # 整合标题和数据
        column_data = zip(columns, [tn, fp, fn, tp, acc, roc, f1, ef, y_pred, y_pred_prob, y_true])
        # 写入CSV
        pd.DataFrame(column_data).T.to_csv(self.model_metric_csv, index=False, mode='a', header=False)

    def plot_roc(self):
        # 读取数据
        df = pd.read_csv(self.model_metric_csv, encoding='utf-8')
        # 获取y_true, y_pred
        y_true, y_pred_proba = eval(df.loc[0, 'y_true']), eval(df.loc[0, 'y_pred_proba'])
        # 获取ROC_AUC
        roc = df.loc[0, 'roc']
        # 获取FPT,RECALL
        FPR, recall, thresholds = roc_curve(y_true, y_pred_proba, pos_label=1)
        # 作图
        plt.figure(figsize=(6, 6))  # 创建画布
        plt.plot(FPR, recall, color='red',
                 label='ROC curve (area = {:0.2f})'.format(roc))  # 画曲线
        plt.plot([0, 1], [0, 1], color='black', linestyle='--')  # 画参考线
        plt.xlim([-0.05, 1.05])  # X坐标范围
        plt.ylim([-0.05, 1.05])  # y坐标范围
        plt.xlabel('False Positive Rate')  # x标签
        plt.ylabel('Recall')  # y标签
        plt.title('Receiver operating characteristic')  # 标题
        plt.legend(loc="lower right")  # 标签
        plt.savefig(self.roc_png, format='png', transparent=True, dpi=600)  # 储存图片
        plt.close()  # 关闭对象，否则shap会和roc图重叠

    def plot_shap(self, train_x, clf, feature_names, model_type):
        # 根据算法类型选择合适的解释器
        if model_type == 'tree':
            # 实例化解释器
            explainer = shap.TreeExplainer(clf)
            # 把特征矩阵转成DataFrame
            train_x = pd.DataFrame(train_x, columns=feature_names)
            # 计算shap值
            shap_values = explainer.shap_values(train_x)
            # 作图
            shap.summary_plot(shap_values, train_x, show=False, plot_size=(12, 6))
            # 布局紧密，保证坐标轴名称能被显示出来
            plt.tight_layout()
            # 保存图
            plt.savefig(self.shap_png, dpi=600, transparent=True)
        elif model_type == 'kernel':  # SVM模型
            pass
            # explainer = shap.KernelExplainer(clf)

    def save_best_parameter(self, fina_parameters):
        # 定义写入配置文件的内容
        config_con = '#best_parameter\n{}\n'.format(fina_parameters)
        # 保存到config中
        self.write_config(config_con)
        # 题目-数据
        parameter_value = [(i, j) for i, j in fina_parameters.items()]
        # 保存到CSV中
        pd.DataFrame(parameter_value).T.to_csv(self.final_parameter_csv, index=False, mode='a', header=False)

    def pred2name(self, df, y_pred_proba, pred_labels=None):
        # 读取名字df
        names = df.iloc[:, 0].tolist()
        # 获取标签为1的概率
        y_pred = [i[1] for i in y_pred_proba]
        # 判断长度是否一致
        assert len(names) == len(y_pred), '长度不一致'
        # 合并df
        df = pd.DataFrame(names, columns=['name'])
        df['score'] = y_pred
        # 根据概率值排序，并且给缺失值填充0
        df.sort_values(by='score', inplace=True, ascending=False)
        df.reset_index(drop=True, inplace=True)
        # 选取概率大于=0.5的小分子
        if pred_labels is None:
            labels = np.zeros(len(df))
            labels[df[df.loc[:, 'score'] >= 0.5].index.tolist()] = 1
        else:
            labels = pred_labels
        df['label'] = labels 
        # 写入CSV
        df.to_csv(self.active_ligand_csv, index=False, mode='a')
 

class xgb_model(ml_model):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 超参数范围
        self.hyper_parameters_range = {
            # 0.1表示取不到0
            'n_estimators': [1, 9999],
            'learning_rate': [0.1, 1],
            'subsample': [0.1, 1],
            'max_depth': [1, 9999],
            'gamma': [0, 9999],
            'min_child_weight': [0, 9999],
            'colsample_bytree': [0.1, 1],
            'colsample_bylevel': [0.1, 1],
            'colsample_bynode': [0.1, 1],
            'reg_alpha': [0, 9999],
            'reg_lambda': [0, 9999]
        }
        # 超参数默认值
        self.hyper_parameters = {
            'n_estimators': 100,
            'learning_rate': 0.3,
            'subsample': 1,
            'max_depth': 6,
            'gamma': 0,
            'min_child_weight': 1,
            'colsample_bytree': 1,
            'colsample_bylevel': 1,
            'colsample_bynode': 1,
            'reg_alpha': 0,
            'reg_lambda': 1
        }
        self.model = XGBClassifier

    def get_xgb(self, best_parameter):  # 获取模型
        return XGBClassifier(n_estimators=int(best_parameter['n_estimators']),
                             learning_rate=float(best_parameter['learning_rate']),
                             subsample=float(best_parameter['subsample']),
                             max_depth=int(best_parameter['max_depth']),
                             gamma=float(best_parameter['gamma']),
                             min_child_weight=float(best_parameter['min_child_weight']),
                             colsample_bytree=float(best_parameter['colsample_bytree']),
                             colsample_bylevel=float(best_parameter['colsample_bylevel']),
                             colsample_bynode=float(best_parameter['colsample_bynode']),
                             reg_alpha=float(best_parameter['reg_alpha']),
                             reg_lambda=float(best_parameter['reg_lambda']),
                             n_jobs=-1,
                             random_state=42)


class rf_model(ml_model):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 超参数范围
        self.hyper_parameters_range = {
            'n_estimators': [1, 9999],
            'max_depth': [1, 9999],
            'min_samples_split': [2, 9999],
            'min_samples_leaf': [1, 9999],
            'max_features': [1, 19999, 'sqrt', 'log2'],
            'min_impurity_decrease': [0, 9999]
        }
        # 超参数默认值
        self.hyper_parameters = {
            'n_estimators': 100,
            'max_depth': 100,
            'min_samples_split': 2,
            'min_samples_leaf': 1,
            'max_features': 'sqrt',
            'min_impurity_decrease': 0
        }
        self.model = RandomForestClassifier

    def get_rf(self, best_parameter):  # 获取模型
        # 对最大特征进行判断
        max_feature = best_parameter['max_features']  # 赋值
        if max_feature.isdigit():  # 若为数字，转换成数字， 否则以字符串传入
            max_feature = int(max_feature)
        return RandomForestClassifier(n_estimators=int(best_parameter['n_estimators']),
                                      max_depth=int(best_parameter['max_depth']),
                                      min_samples_leaf=int(best_parameter['min_samples_leaf']),
                                      min_samples_split=int(best_parameter['min_samples_split']),
                                      max_features=max_feature,
                                      min_impurity_decrease=float(best_parameter['min_impurity_decrease']),
                                      n_jobs=-1,
                                      random_state=42)


class svm_model(ml_model):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 超参数范围
        self.hyper_parameters_range = {
            # 0.1表示取不到0
            'C': [0.1, 9999],
            'gamma': [0.1, 9999]
        }
        # 超参数默认值
        self.hyper_parameters = {
            'C': 1,
            'gamma': 'auto'
        }
        self.model = SVC

    def get_svm(self, best_parameter):  # 获取模型
        # 对gamma进行格式转换
        gamma = best_parameter['gamma']
        # 判断是否为数字
        try:
            gamma = float(gamma)  # 若为数字，转换成浮点数
        except:
            gamma = 'auto'  # 若不为数字，恢复默认参数
        # 返回
        return SVC(C=float(best_parameter['C']),
                   gamma=gamma,
                   class_weight='balanced',
                   random_state=42, probability=True)
