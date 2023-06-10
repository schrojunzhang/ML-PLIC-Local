#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-10

import os
from base_scripts.base_class.soft_base import soft

class interaction_fingerprint(soft):

    def __init__(self, job_id):
        super().__init__(job_id)
        self.ifp2fingerprint = {
            'ifp': ['name'] + ['ifp_{}'.format(i) for i in range(1, 1025)],
            'splif': ['name'] + ['splif_{}'.format(i) for i in range(1, 1025)],
            'plec': ['name'] + ['plec_{}'.format(i) for i in range(1, 1025)],
            'silirid': ['name'] + ['silirid_{}'.format(i) for i in range(1, 169)]
        }

    def cal_ifp(self, ligand):
        # 指纹类型
        fp_type = 'ifp'
        # 还原路径
        path_files = self.path_for_lig  # 分子存放路径
        path_results = self.path_for_csv  # csv存放路径
        # 蛋白文件
        protein_file = self.protein
        # 定义小分子
        ligand = '{}/{}'.format(path_files, ligand)
        # 结果CSV文件
        csv_file = '{}/{}.csv'.format(path_results, fp_type)  # csv文件
        # 计算
        self.oddt_excute(protein_file, ligand, csv_file, fp_type)

    def cal_splif(self, ligand):
        # 指纹类型
        fp_type = 'splif'
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        path_results = self.path_for_csv  # csv存放路径
        # 蛋白文件
        protein_file = self.protein
        # 定义小分子
        ligand = '{}/{}'.format(path_files, ligand)
        # 结果CSV文件
        csv_file = '{}/{}.csv'.format(path_results, fp_type)  # csv文件
        # 计算
        self.oddt_excute(protein_file, ligand, csv_file, fp_type)

    def cal_plec(self, ligand):
        # 指纹类型
        fp_type = 'plec'
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        path_results = self.path_for_csv  # csv存放路径
        # 蛋白文件
        protein_file = self.protein
        # 定义小分子
        ligand = '{}/{}'.format(path_files, ligand)
        # 结果CSV文件
        csv_file = '{}/{}.csv'.format(path_results, fp_type)  # csv文件
        # 计算
        self.oddt_excute(protein_file, ligand, csv_file, fp_type)

    def cal_silirid(self, ligand):
        # 指纹类型
        fp_type = 'silirid'
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        path_results = self.path_for_csv  # csv存放路径
        # 蛋白文件
        protein_file = self.protein
        # 定义小分子
        ligand = '{}/{}'.format(path_files, ligand)
        # 结果CSV文件
        csv_file = '{}/{}.csv'.format(path_results, fp_type)  # csv文件
        # 计算
        self.oddt_excute(protein_file, ligand, csv_file, fp_type)

class fingerpringt(soft):

    def __init__(self, job_id):
        super().__init__(job_id)
        self.fp2fingerprint = {
            'ecfp': ['name'] + ['ecfp_{}'.format(i) for i in range(1, 1025)],
            'pubchem': ['name'] + ['pubchem_{}'.format(i) for i in range(1, 882)],
        }

    def cal_ecfp(self, ligand):
        # 指纹类型
        fp_type = 'ecfp'
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        path_results = self.path_for_csv  # csv存放路径
        dst_path = '{}/ecfp'.format(path_job)  # 对应路径
        # 定义小分子
        lig_name = ligand.split('.')[0]  # 分子纯名字
        ligand_file = '{}/{}'.format(path_files, ligand)  # 分子全路径
        pre_ligand = '{}/{}.mol'.format(dst_path, lig_name)  # 准备后的分子
        # 准备分子,转换格式
        self.openeye_convert(ligand_file, pre_ligand)
        # 结果CSV文件
        csv_file = '{}/{}.csv'.format(path_results, fp_type)  # csv文件
        # 计算
        self.rdkit_excute(pre_ligand, csv_file, fp_type)

    def cal_pubchem(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/pubchem'.format(path_job)  # 指纹文件夹
        # 定义小分子
        lig_name = ligand.split('.')[0]  # 小分子纯名
        ligand_file = '{}/{}'.format(path_files, ligand)  # 分子全路径
        # 日志文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 计算
        cmd = 'cd /home/xujun/Soft/Score_Function/PaDel &&java -jar PaDEL-Descriptor.jar -dir {} -waitingjobs ' \
              '-1 -fingerprints -file {}'.format(ligand_file, log_file)
        # 执行
        os.system(cmd)


