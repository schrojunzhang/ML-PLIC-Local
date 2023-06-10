#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-15


import sys
import pandas as pd
from multiprocessing import Pool
from base_scripts.base_class.docking_base import vina_docking

# 获取任务名
job_id = sys.argv[1]
# 实例化对象
vina = vina_docking(job_id)
# 获取表单
job_form = vina.read_config(job_type=vina.job_type)
# 从表单中组装出命令行
protein_cmd, ligand_cmd, docking_cmd = vina.form2cmd(job_form)
# 蛋白准备
vina.protein_prep(protein_cmd)
# 获取小分子
ligands = vina.ligands
# 创建CSV
pd.DataFrame(['name', 'score']).T.to_csv(vina.docking_score_csv, index=False, header=False)
# 组成列表
cmd_ligands = zip([ligand_cmd for i in range(len(ligands))], ligands)  # 小分子列表
# 多进程准备并对接小分子
pool = Pool(28)
# 准备小分子
print('start ligand preparation')
pool.starmap(vina.ligand_prep, cmd_ligands)
# 判断蛋白活性位点
print('get binding pocket')
x, y, z = vina.read_pocket()
# 判断使用的对接软件
soft = job_form[vina.soft_form]
print('using docking program: {}'.format(soft))
# 组成列表
cmd_soft_ligands_pdbqt = zip([docking_cmd % (x, y, z) for i in range(len(ligands))], [soft for i in range(len(ligands))],
                        ligands)  # 小分子列表
# 对接
print('start docking')
pool.starmap(vina.ligand_docking, cmd_soft_ligands_pdbqt)  # 小分子对接打分
# 搜集数据写入CSV
print('write score to CSV file')
pool.map(vina.get_score, ligands)
pool.close()
pool.join()
# 转换蛋白格式到描述符提取文件夹
vina.openbabel_transform(src_file=vina.protein_pdbqt, dst_file=vina.protein_docked)
