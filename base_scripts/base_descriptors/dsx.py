#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-09

import os
import glob
import sys
import shutil
import pandas as pd
from base_scripts.base_class import scoring_function_base
from multiprocessing import Pool

def main():
    sf = sys.argv[0].split('/')[-1].split('.')[0]
    job_id = sys.argv[1]
    # 实例化打分函数
    sf_py = scoring_function_base.scoring_function(job_id)
    # 还原路径
    path_job = sf_py.path_job  # 任务文件夹路径
    path_results = sf_py.path_for_csv  # csv存放路径
    dst_path = '{}/{}'.format(path_job, sf)  # 目标路径
    # 获取分子
    ligands = sf_py.ligands  # 小分子名字合集
    # 多进程计算
    pool = Pool(28)
    pool.map(sf_py.cal_dsx, ligands)
    pool.close()
    pool.join()
    # 收集数据
    energy_term = sf_py.sf2energy_term[sf]  # 能量项
    csv_file = '{}/{}.csv'.format(path_results, sf)  # csv文件
    pd.DataFrame(energy_term).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 循环读取日志获取数据，添加至CSV中
    for ligand in ligands:
        lig_name = ligand.split('.')[0]  # 获取纯名字
        result = [lig_name]  # 初始化空列表
        for i in [0, 2, 3]:  # 循环读取三种能量项txt
            tag = len(result)
            log = '{}/{}_{}.txt'.format(dst_path, lig_name, i)  # 定义日志文件
            if os.path.exists(log):  # 若文件存在，则读取数据
                # 读取数据
                with open(log, 'r') as f:
                    con = f.readlines()
                # 整理数据
                for j in range(len(con)):  # 循环读取行
                    try:  # 预防下标越界错误
                        if con[j].startswith('@RESULTS'):
                            result.append(float(con[j + 4].split('|')[3].strip()))
                            break
                    except:
                        pass
                # 防止错位
                if len(result) != tag + 1:
                    result.append(' ')
        # 写数据到CSV中
        pd.DataFrame(result).T.to_csv(csv_file, index=False, header=False, mode='a')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()
