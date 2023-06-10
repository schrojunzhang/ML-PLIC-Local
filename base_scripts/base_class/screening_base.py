#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-17

from base_scripts.base_class.global_base import screening_control
from base_scripts.base_class.model_base import svm_model, xgb_model, rf_model


class screening_workflow(screening_control, svm_model, xgb_model, rf_model):

    def __init__(self, job_id):
        super().__init__(job_id)
        super(xgb_model, self).__init__(job_id)
        super(svm_model, self).__init__(job_id)
        super(rf_model, self).__init__(job_id)
        # 获取表单
        # 表单
        self.dock_form = self.read_config(job_type=self.dock_job_type)
        self.descirptor_form = self.read_config(job_type=self.descriptor_job_type)
        self.model_form = self.read_config(job_type=self.model_job_type)
        self.best_parameters = self.read_config(job_type=self.para_prefix)
        # 从表单获取模型对象
        self.form2model = {
            'svm': svm_model(job_id),
            'xgb': xgb_model(job_id),
            'rf': rf_model(job_id),
        }
        # 从表单获取分类器对象
        self.form2clf = {
            'svm': svm_model(job_id).get_svm,
            'xgb': xgb_model(job_id).get_xgb,
            'rf': rf_model(job_id).get_rf
        }
