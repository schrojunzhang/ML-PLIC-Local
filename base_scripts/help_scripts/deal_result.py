#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-18
# 用于数据结果的处理
import sys
from base_scripts.base_class.queue_result_base import result
from base_scripts.base_class.global_base import pipline_control

def main(job_id, job_type):
    # 实例化对象
    this_result = result(job_id)
    # 根据不同任务类型做不同处理
    # 闭包函数
    def final_result(job_type):
        # 获取得分前100的小分子，组成复合物
        active_names = this_result.top100_ligs(this_result.type2ligandcsv[job_type])  # 获取名字
        this_result.active2complex(active_names)  # 组成复合物
        # 计算聚类
        print('start clustering')
        this_result.cal_cluster(active_names)
        # 计算相互作用
        print('start calculating interaction')
        this_result.cal_interactions(job_id)
        # 如果是对接任务，合并分子
        if job_type == 'docking':
            this_result.merge_lig()
    # 对接处理
    if job_type == 'docking':
        final_result(job_type)
        # # 获取得分前100的小分子，组成复合物
        # active_names = this_result.top100_ligs(this_result.type2ligandcsv[job_type])  # 获取名字
        # this_result.active2complex(active_names)  # 组成复合物
        # # 计算聚类
        # this_result.cal_cluster(active_names)
        # # 计算相互作用
        # this_result.cal_interactions(job_id)
    elif job_type == 'pipeline':
        # 如果是pipeline， 进一步判断是否进行对接和建模
        this_pipeline = pipline_control(job_id)
        docking = this_pipeline.whether_dock()
        modelling = this_pipeline.whether_modelling()
        # 如果不对接也不建模
        if not modelling and not docking:
            pass
        else:
            # 对接，建模至少其一
            if not modelling:
                # 若只对接
                job_type = 'docking'
            final_result(job_type)
            # # 获取得分前100的小分子，组成复合物
            # active_names = this_result.top100_ligs(this_result.type2ligandcsv[job_type])  # 获取名字
            # this_result.active2complex(active_names)  # 组成复合物
            # # 计算聚类
            # this_result.cal_cluster(active_names)
            # # 计算相互作用
            # this_result.cal_interactions(job_id)
    # 压缩文件
    print('start zipping')
    this_result.zip_result()


if __name__ == '__main__':
    # 获取任务名
    job_id = sys.argv[1]
    # 任务类型
    job_type = sys.argv[2]
    # 执行
    main(job_id, job_type)
