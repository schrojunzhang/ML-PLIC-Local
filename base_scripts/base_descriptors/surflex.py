#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-10
# 先把蛋白.pdb 转成 mol2， 然后为每个分子创建文件夹， 通过surflex_getpocket获取口袋， 生成p1-protomol.mol2
# 创建list文件，内容为小分子文件路径， 通过指令 输入 ./list ./p1-protomol.mol2 蛋白文件路径  日志路径
# 生成了如下文件： scores， loggeom（日志文件）, loggeom-results.mol2

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
    # 预处理蛋白
    pre_protein = '{}/protein_file.mol2'.format(dst_path)
    sf_py.openeye_convert(protein, pre_protein)
    # 共晶配体
    crystal_ligand = sf_py.crystal_ligand
    # 获取口袋生成p1-protomol.mol2文件
    sf_py.surflex_getpocket(crystal_ligand, pre_protein, dst_path)
    # 多进程计算
    pool = Pool(28)
    pool.map(sf_py.cal_surflex, ligands)
    pool.close()
    pool.join()
    # 收集数据
    energy_term = sf_py.sf2energy_term[sf]  # 能量项
    csv_file = '{}/{}.csv'.format(path_results, sf)  # csv文件
    pd.DataFrame(energy_term).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 循环读取日志获取数据，添加至CSV中
    for ligand in ligands:
        lig_name = ligand.split('.')[0]  # 获取纯名字
        log = '{}/{}/loggeom'.format(dst_path, lig_name)
        result = []  #初始化空列表
        if os.path.exists(log):  # 若文件存在，则读取数据
            # 读取数据
            with open(log, 'r') as f:
                con = f.readlines()
            # 整理数据
            tem_list = con[0].split()  # rot + total + crash + self_crash + polar
            result =[tem_list[3].strip()] + [tem_list[13+i].strip() for i in range(1, 7, 2)] + [tem_list[20].rstrip(']').strip()]
        result.insert(0, lig_name)  # 插入分子名称
        # 写数据到CSV中
        pd.DataFrame(result).T.to_csv(csv_file, index=False, header=False, mode='a')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()