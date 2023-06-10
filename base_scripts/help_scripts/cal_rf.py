#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-09

import os
import sys
import glob
import csv
from rfscore.ob import get_molecule
from rfscore.credo import contacts
from multiprocessing import Pool
from multiprocessing import Manager


def rfscore_calc(protein, ligand, mydesc):
    if mydesc == 'element':
        descriptor, labels = contacts.element_descriptor(protein, ligand, binsize=1)
    elif mydesc == 'sybyl':
        descriptor, labels = contacts.sybyl_atom_type_descriptor(protein, ligand, binsize=1)
    elif mydesc == 'credo':
        descriptor, labels = contacts.sift_descriptor(protein, ligand, binsize=1)
    return descriptor, labels


def rfscore_calc_parallel(lig_name, ligand_file, score_type, i, return_dict):
    protein_molecule = get_molecule(protein)
    ligand_molecule = get_molecule(ligand_file)
    descriptor, labels = rfscore_calc(protein_molecule, ligand_molecule, score_type)
    descriptor_list = descriptor.tolist()
    descriptor_list.insert(0, lig_name)
    return_dict[i] = descriptor_list

# 任务文件夹路径
path_job = sys.argv[1]
# 目标路径
dst_path = sys.argv[2]
# 打分类型
score_type = sys.argv[3]
# 分子存放路径
path_files = '{}/files'.format(path_job)
# csv存放路径
path_results = '{}/csvs'.format(path_job)
# 获取分子
ligands = [i.split('/')[-1] for i in glob.glob('{}/*.mol2'.format(path_files))]
protein = '{}/protein_file.pdb'.format(path_files)
# 结果文件
csv_file = '{}/rfscore_{}.csv'.format(path_results, score_type)
# 多线程计算
manger = Manager()  # 实例化对象
return_dict = manger.dict()
pool = Pool(28)
jobs = []
for ligand in ligands:
    j = ligands.index(ligand)
    ligand_file = '{}/{}'.format(path_files, ligand)  # 小分子文件
    lig_name = ligand.split('.')[0]  # 小分子纯名字
    p = pool.apply_async(rfscore_calc_parallel, (lig_name, ligand_file, score_type, j, return_dict))
    jobs.append(p)
pool.close()
pool.join()
# 整合结果到CSV中
results = return_dict.values()  # 结果列表
mycsv = open(csv_file, 'a')  # 打开CSV
mycsvwriter = csv.writer(mycsv)  # 初始化对象
for result in results:  # 循环写入数据
    mycsvwriter.writerow(result)
mycsv.close()