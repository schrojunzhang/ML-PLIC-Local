#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-10

import os
import glob
import sys
import shutil
import numpy as np
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
    ligands = sf_py.ligands
    protein = sf_py.protein
    # 获取坐标
    x, y, z = sf_py.read_pocket()
    # 预处理蛋白质
    protein_temp_file = '{}/protein_temp.pdb'.format(dst_path)
    pre_protein_file = '{}/protein_file.pdbqt'.format(dst_path)
    # 处理蛋白
    if not os.path.exists(pre_protein_file):
        cmdline = 'cat %s | sed \'/HETATM/\'d > %s &&' % (protein, protein_temp_file)
        cmdline += 'module purge &&'
        cmdline += 'module load autodock &&'
        cmdline += 'prepare_receptor4.py -r %s -o %s -A hydrogens -U nphs_lps_waters_nonstdres &&' % (
            protein_temp_file, pre_protein_file)
        cmdline += 'rm  %s' % protein_temp_file
        os.system(cmdline)
    # 列表
    length = len(ligands)
    ligands_xyz = zip(ligands, [x] * length, [y] * length, [z] * length)
    # 多进程计算
    pool = Pool(28)
    pool.starmap(sf_py.cal_vina, ligands_xyz)
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
            for i in range(len(con)):
                if con[i].startswith('Affinity:') and len(con) - i >= 7:
                    score = con[i].strip().split(':')[1].split('(')[0].strip()
                    gauss1 = con[i + 2].split(':')[-1].strip()
                    gauss2 = con[i + 3].split(':')[-1].strip()
                    repulsion = con[i + 4].split(':')[-1].strip()
                    hydrophobic = con[i + 5].split(':')[-1].strip()
                    hydrogen_bond = con[i + 6].split(':')[-1].strip()
                    if float(score) != 0:
                        rt = (((-0.035579) * float(gauss1) + (-0.005156) * float(gauss2) + (
                            0.84024500000000002) * float(
                            repulsion) + (-0.035069000000000003) * float(hydrophobic) + (-0.58743900000000004) * float(
                            hydrogen_bond)) / float(score) - 1) // 0.058459999999999998
                    else:
                        rt = np.nan
                    result = [score, gauss1, gauss2, repulsion, hydrophobic, hydrogen_bond, rt]
                    break
        result.insert(0, lig_name)  # 插入分子名称
        # 写数据到CSV中
        pd.DataFrame(result).T.to_csv(csv_file, index=False, header=False, mode='a')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()