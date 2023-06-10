#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-05

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
    # 转换蛋白格式
    pre_protein_file = '{}/protein_file.mae'.format(dst_path)
    cmd = 'cd {}&&module load schrodinger&&structconvert {} {}'.format(dst_path, protein, pre_protein_file)
    os.system(cmd)
    # 准备生成格点
    # 蛋白加电荷
    sf_py.dock_pro_pre(pre_protein_file, dst_path)
    # 获取活性口袋
    sf_py.dock_active_site(dst_path, crystal_ligand)
    # 生成格点
    sf_py.dock_generate_grid(dst_path)
    # 多进程计算
    pool = Pool(28)
    pool.map(sf_py.cal_hawkins, ligands)
    pool.close()
    pool.join()
    # 收集数据
    energy_term = sf_py.sf2energy_term[sf]  # 能量项
    csv_file = '{}/{}.csv'.format(path_results, sf)  # csv文件
    pd.DataFrame(energy_term).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 循环读取日志获取数据，添加至CSV中
    for ligand in ligands:
        lig_name = ligand.split('.')[0]  # 获取纯名字
        log = '{}/{}_scored.mol2'.format(dst_path, lig_name)  # 定义日志文件
        result = [lig_name]  #初始化空列表
        if os.path.exists(log):  # 若文件存在，则读取数据
            # 读取数据
            with open(log, 'r') as f:
                con = f.readlines()
            # 整理数据
            for i in range(len(con)):
                try:  # 预防下标越界错误
                    if con[i].startswith('##########                                Name:'):
                        result = [lig_name] + [con[i+j].split(':')[-1].strip() for j in range(1, 7)]
                        break
                except:
                    pass
        # 写数据到CSV中
        pd.DataFrame(result).T.to_csv(csv_file, index=False, header=False, mode='a')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()
