#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-10

import sys
import shutil
import glob
import pandas as pd
from base_scripts.base_class import fingerprint_base
from multiprocessing import Pool

def main():
    fp_type = sys.argv[0].split('/')[-1].split('.')[0]
    job_id = sys.argv[1]
    # 实例化打分函数
    ifp = fingerprint_base.interaction_fingerprint(job_id)
    # 还原路径
    path_job = ifp.path_job  # 任务文件夹路径
    path_files = ifp.path_for_lig  # 分子存放路径
    path_results = ifp.path_for_csv  # csv存放路径
    dst_path = '{}/{}'.format(path_job, fp_type)  # 目标路径
    # 获取分子
    ligands = ifp.ligands
    # csv文件标题
    header = ifp.ifp2fingerprint[fp_type]  # 标题
    csv_file = '{}/{}.csv'.format(path_results, fp_type)  # csv文件
    pd.DataFrame(header).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 多进程计算
    pool = Pool(28)
    pool.map(ifp.cal_ifp, ligands)
    pool.close()
    pool.join()
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()