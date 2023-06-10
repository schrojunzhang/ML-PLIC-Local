#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-15

import os
import sys
from base_scripts.base_class.model_base import rf_model


def main():
    # 获取任务名
    job_id = sys.argv[1]
    # 实例化对象
    this_model = rf_model(job_id)
    # 获取表单
    job_form = this_model.read_config(job_type=this_model.job_type)
    # 获取df
    df = this_model.get_df()
    # 舍弃空值
    df.dropna(inplace=True)
    # 获取x,y 获取特征名字
    x, y, feature_names = this_model.get_x_y(df)
    # 划分训练集和测试集
    train_size = float(job_form['train_size'])  # 获取比例
    train_x, test_x, train_y, test_y = this_model.train_test_split(x=x, y=y, train_size=train_size)  # 划分
    # 数据预处理  #######################################################################################################
    # 缩放
    scaler = job_form['scaler']  # 从表单获取缩放对象
    train_x, test_x = this_model.data_scaler(scaler, train_x=train_x, test_x=test_x)  # 数据缩放
    # 低方差过滤
    train_x, test_x, feature_names = this_model.remove_same_feature(train_x=train_x, test_x=test_x, feature_names=feature_names)
    # 特征选择
    selector = job_form['feature_selector']  # 获取特征选择器
    selector_para = job_form['feature_selector_para']  # 选择器参数
    train_x, test_x, feature_names = this_model.feature_selection(selector, train_x=train_x, train_y=train_y, test_x=test_x,
                                                   para=selector_para, feature_names=feature_names)  # 特征选择
    # 输出特征重要性
    this_model.cal_feature_importance(feature_names, x=train_x, y=train_y)
    # 获得寻优次数  #####################################################################################################
    max_eval = int(job_form['max_eval'])
    # 获取超参数
    parameter_space = this_model.get_space(job_form=job_form, hyper_parameters=this_model.hyper_parameters,
                                           hyper_parameters_range=this_model.hyper_parameters_range)
    # 判断是否进行寻优
    if max_eval > 0 and parameter_space:  # 进行
        # 获得模型对象，用于超参数寻优
        model = this_model.model
        # 寻优
        best_parameters, trials = this_model.hyperopt_tuning(train_x=train_x, train_y=train_y, model=model, max_evals=max_eval,
                                                     parameter_space=parameter_space)
        # 学习曲线图
        this_model.plot_learning_curve()
        # 方差偏差分析 作图
        this_model.bias_var_analysis()
        # 超参数寻优图
        this_model.plot_hyperParameter_scatter(trials)
        # 根据最优超参数及参数空间组合组成分类器
        fina_parameters = this_model.best2parameter(job_form=job_form, hyper_parameters=this_model.hyper_parameters,
                                                    best_parameters=best_parameters)
    else:  # 不进行寻优
        fina_parameters = this_model.hyper_parameters
    #  获取最优模型  ####################################################################################################
    clf = this_model.get_rf(fina_parameters)
    # 训练
    clf.fit(train_x, train_y)
    # 预测
    y_pred = clf.predict(test_x)
    # 预测概率
    y_pred_proba = clf.predict_proba(test_x)
    # 模型评价
    this_model.model_metric(y_true=test_y, y_pred=y_pred, y_pred_prob=y_pred_proba)
    # 绘制ROC曲线图
    this_model.plot_roc()
    # 保存参数
    this_model.save_best_parameter(fina_parameters)
    try:
        # 绘制shap模型解释图
        this_model.plot_shap(train_x, clf, feature_names=feature_names, model_type='tree')
    except:
        pass

if __name__ == '__main__':
    main()

