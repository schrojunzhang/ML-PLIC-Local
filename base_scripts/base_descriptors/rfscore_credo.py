#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-09

import sys
import shutil
import pandas as pd
from base_scripts.base_class import scoring_function_base


def main():
    sf = sys.argv[0].split('/')[-1].split('.')[0]
    job_id = sys.argv[1]
    # 实例化打分函数
    sf_py = scoring_function_base.scoring_function(job_id)
    # 还原路径
    path_job = sf_py.path_job  # 任务文件夹路径
    path_results = sf_py.path_for_csv  # csv存放路径
    dst_path = '{}/{}'.format(path_job, sf)  # 目标路径
    # csv文件标题
    energy_term = sf_py.sf2energy_term[sf]
    csv_file = '{}/{}.csv'.format(path_results, sf)  # csv文件
    pd.DataFrame(energy_term).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 计算
    sf_py.rfscore_excute(path_job, dst_path, 'credo')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()