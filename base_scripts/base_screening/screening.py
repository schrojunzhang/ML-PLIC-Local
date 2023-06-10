#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-16

import sys
from base_scripts.base_class.screening_base import screening_workflow


def main():
    # 获取ID
    job_id = sys.argv[1]
    # 实例化对象
    this_job = screening_workflow(job_id)
    # 获取表单
    job_form = this_job.model_form
    # 根据表单判断导入什么算法模型
    algo = this_job.algorithm[job_form.get('algorithm')]  # 从表单获取算法
    this_model = this_job.form2model[algo]  # 根据算法选取模型
    # 获取df
    train_df = this_model.get_df(csv_file=this_job.descriptor_csv)  # 获取训练集数据
    test_df = this_model.get_df(csv_file=this_job.test_descriptors_csv)  # 获取测试集数据
    # 舍弃空值
    train_df.dropna(inplace=True)
    test_df.dropna(inplace=True)
    # 获取x,y
    train_x, train_y, feature_names = this_model.get_x_y(train_df)  # 划分xy
    test_x, test_y, feature_names = this_model.get_x_y(test_df)
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
    # 获取最佳参数
    fina_parameters = this_job.best_parameters
    #  获取最优模型  ####################################################################################################
    clf = this_job.form2clf[algo](fina_parameters)
    # 训练
    clf.fit(train_x, train_y)
    # 预测
    y_pred_proba = clf.predict_proba(test_x)
    y_pred = clf.predict(test_x)
    # 转换为活性分子名并写入CSV  
    this_job.pred2name(test_df, y_pred_proba, y_pred)


if __name__ == '__main__':
    main()
