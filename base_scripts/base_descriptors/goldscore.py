#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-08

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
    # 获取坐标
    x, y, z = sf_py.read_pocket()
    # 获取分子
    ligands = sf_py.ligands  # 小分子名字合集
    # 列表
    length = len(ligands)
    ligands_xyz = zip(ligands, [x]*length, [y]*length, [z]*length)
    # 多进程计算
    pool = Pool(28)
    pool.starmap(sf_py.cal_goldscore, ligands_xyz)
    pool.close()
    pool.join()
    # 收集数据
    energy_term = sf_py.sf2energy_term[sf]  # 能量项
    csv_file = '{}/{}.csv'.format(path_results, sf)  # csv文件
    pd.DataFrame(energy_term).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 循环读取日志获取数据，添加至CSV中
    for ligand in ligands:
        lig_name = ligand.split('.')[0]  # 配体无后缀名字
        log = '{}/{}/log.csv'.format(dst_path, lig_name)  # 日志文件
        result = []  # 初始化结果列表
        if os.path.exists(log):  # 若存在文件，读取数据
            # 读取数据
            with open(log, 'r') as f:
                con = f.readlines()
            # 整理数据
            try:  # 预防错误
                result = con[1].split(',')[4:10]
            finally:
                result.insert(0, lig_name)
        else:  # 若不存在log文件，则只输出分子名字
            result.append(lig_name)
        # 写数据到CSV中
        pd.DataFrame(result).T.to_csv(csv_file, index=False, header=False, mode='a')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()