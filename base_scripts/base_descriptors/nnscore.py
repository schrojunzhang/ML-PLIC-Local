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
    path_files = sf_py.path_for_lig  # 分子存放路径
    path_results = sf_py.path_for_csv  # csv存放路径
    dst_path = '{}/{}'.format(path_job, sf)  # 目标路径
    # 获取分子
    ligands = sf_py.ligands
    protein = sf_py.protein
    # 预处理蛋白质
    protein_temp_file = '{}/protein_temp.pdb'.format(dst_path)
    pre_protein_file = '{}/protein_file.pdbqt'.format(dst_path)
    if not os.path.exists(pre_protein_file):
        cmd = 'cat {0} | sed \'/HETATM/\'d > {1} &&module purge&&module load vina&&prepare_receptor4.py -r {1} -o {2} ' \
              '-A hydrogens -U nphs_lps_waters_nonstdres &&rm {1}'.format(protein, protein_temp_file, pre_protein_file)
        os.system(cmd)
    # 多进程计算
    pool = Pool(28)
    pool.map(sf_py.cal_nnscore, ligands)
    pool.close()
    pool.join()
    # 收集数据
    # csv文件标题
    energy_term = ['name'] + sf_py.sf2energy_term['nn_vina'] \
                  + ['atp2_%s' % i for i in sf_py.sf2energy_term['nn_close']] \
                  + ['atp4_%s' % i for i in sf_py.sf2energy_term['nn_semi']] \
                  + ['lat_%s' % i for i in sf_py.sf2energy_term['nn_atm']] \
                  + ['ele_%s' % i for i in sf_py.sf2energy_term['nn_elec']] \
                  + sf_py.sf2energy_term['nn_rot'] \
                  + ['siteflex_%s' % i for i in sf_py.sf2energy_term['nn_active_site_flexibility']] \
                  + ['hbond_%s' % i for i in sf_py.sf2energy_term['nn_hb']] \
                  + ['hydrophobic_%s' % i for i in sf_py.sf2energy_term['nn_hydrophobics']] \
                  + ['stacking_%s' % i for i in sf_py.sf2energy_term['nn_stacking']] \
                  + ['pi_cation_%s' % i for i in sf_py.sf2energy_term['nn_pi_cation']] \
                  + ['t_shaped_%s' % i for i in sf_py.sf2energy_term['nn_t_shape']] \
                  + ['salt_bridges_%s' % i for i in sf_py.sf2energy_term['nn_salt_bridges']]
    csv_file = '{}/{}.csv'.format(path_results, sf)  # csv文件
    pd.DataFrame(energy_term).T.to_csv(csv_file, index=False, header=False)  # 首次创建CSV文件
    # 循环读取日志获取数据，添加至CSV中
    for ligand in ligands:
        lig_name = ligand.split('.')[0]  # 获取纯名字
        log = '{}/{}.txt'.format(dst_path, lig_name)  # 定义日志文件
        result = [lig_name]  #初始化空列表
        if os.path.exists(log):  # 若文件存在，则读取数据
            # 读取数据
            with open(log, 'r') as f:
                con = f.readlines()
            # 整理数据
            result += eval(con[0].strip())
        # 写数据到CSV中
        pd.DataFrame(result).T.to_csv(csv_file, index=False, header=False, mode='a')
    # 删除文件夹
    shutil.rmtree(dst_path)

if __name__ == '__main__':
    main()