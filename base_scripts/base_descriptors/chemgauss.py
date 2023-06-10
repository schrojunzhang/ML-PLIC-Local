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
    # 获取分子
    ligands = sf_py.ligands  # 小分子名字合集
    protein = sf_py.protein  # 蛋白文件
    crystal_ligand = sf_py.crystal_ligand  # 晶体结构
    # 预处理蛋白质
    protein_temp_file = '{}/protein_file.mol2'.format(dst_path)
    pre_protein_file = '{}/protein_file.oeb'.format(dst_path)
    # 处理蛋白
    cmdline = 'cd {} &&'.format(dst_path)
    cmdline += 'module load openeye &&'
    cmdline += 'convert.py {} {} &&'.format(protein, protein_temp_file)
    cmdline += 'export PATH=$PATH:/opt/openeye/applications-2018.11.3/bin &&'
    cmdline += 'receptor_setup -protein {} -bound_ligand {} -receptor {}&&'.format(protein_temp_file, crystal_ligand, pre_protein_file)
    cmdline += 'rm {}'.format(protein_temp_file)
    os.system(cmdline)
    # 多进程计算
    pool = Pool(28)
    pool.map(sf_py.cal_chemgauss, ligands)
    pool.close()
    pool.join()
    # 收集数据
    energy_term = sf_py.sf2energy_term[sf]  # 能量项
    csv_file = '{}/{}.csv'.format(path_results, sf)  # csv文件
    pd.DataFrame(energy_term).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 循环读取日志获取数据，添加至CSV中
    for ligand in ligands:
        lig_name = ligand.split('.')[0]  # 获取纯名字
        log = '{}/{}.txt'.format(dst_path, lig_name)  # 定义日志文件
        result = []  #初始化空列表
        if os.path.exists(log):  # 若文件存在，则读取数据
            # 读取数据
            with open(log, 'r') as f:
                con = f.readlines()
            # 整理数据
            if len(con) == 2:
                result = [con[1].split()[i].strip() for i in range(1, 8)]
        result.insert(0, lig_name)  # 插入分子名称
        # 写数据到CSV中
        pd.DataFrame(result).T.to_csv(csv_file, index=False, header=False, mode='a')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()