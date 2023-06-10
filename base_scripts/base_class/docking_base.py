#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-15

import os
import shutil
import pandas as pd
import glob
from base_scripts.base_class.global_base import docking_control


class vina_docking(docking_control):

    def __init__(self, job_id):
        super().__init__(job_id=job_id)
        # 表单组装命令行
        self.soft_form = 'whether_dock'
        self.protein_form = ['protein_repair', 'protein_delchain', 'protein_cleanup']
        self.ligand_form = ['ligand_repair', 'ligand_charge', 'ligand_cleanup']
        self.docking_form = ['space_size', 'exhaustiveness', 'num_modes']
        self.form2para = {
            'protein_repair': {'bonds_hydrogens': ' -A bonds_hydrogens',
                               'bonds': ' -A bonds',
                               'hydrogens': ' -A hydrogens',
                               'checkhydrogens': ' -A checkhydrogens',
                               'None': ''},
            'protein_charge': {'yes': '',
                               'no': ' -C'},
            'protein_delchain': {'yes': ' -e',
                                 'no': ''},
            'protein_cleanup': {'standard': ' -U nphs_lps_waters_nonstdres',
                                'non-polar hydrogens': ' - U nphs',
                                'lone pairs': ' -U lps',
                                'water residues': ' -U waters',
                                'chains': ' -U nonstdres'},
            'ligand_repair': {'bonds_hydrogens': ' -A bonds_hydrogens',
                              'bonds': ' -A bonds',
                              'hydrogens': ' -A hydrogens',
                              'None': ''},
            'ligand_charge': {'yes': '',
                              'no': ' -C'},
            'ligand_cleanup': {'standard': ' -U nphs_lps',
                               'non-polar hydrogens': ' -U nphs',
                               'lone pairs': ' -U lps'},
            'space_size': ' --size_x {0} --size_y {0} --size_z {0}',
            'exhaustiveness': ' --exhaustiveness {}',
            'num_modes': ' --num_modes {}'  
        }
        self.soft2cmd = {'QuickVina': '/home/xujun/Soft/Score_Function/qvina-master/bin/qvina2.1',
                         'Smina': '/home/xujun/Soft/Score_Function/smina/smina --addH False -q'}
        # 文件
        self.docking_log = self.path_for_dock + '/{}.log'  # 单个分子的打分log

    def form2cmd(self, job_form):  # 从表单组装命令行
        # 初始化字符串
        protein_prep = 'prepare_receptor4.py -r {} -o {}'.format(self.protein, self.protein_pdbqt)
        ligand_prep = 'prepare_ligand4.py -l {} -o {}'
        docking_cmd = self.soft2cmd[
                          job_form[self.soft_form]] + ' --receptor {} --ligand {} --out {} --log {} --center_x %s ' \
                                                      '--center_y %s --center_z %s --cpu 28'
        # 从表单获取参数
        # 蛋白质
        for item in self.protein_form:
            protein_prep += self.form2para[item][job_form[item]]
        # 小分子
        for item in self.ligand_form:
            ligand_prep += self.form2para[item][job_form[item]]
        # 对接
        for item in self.docking_form:
            # 根据参数不同，转换数字格式
            if item == 'space_size':
                parameter = float(job_form[item])
            else:
                parameter = int(job_form[item])
            # 添加参数
            docking_cmd += self.form2para[item].format(parameter)
        # 返回
        return protein_prep, ligand_prep, docking_cmd

    def protein_prep(self, cmd):
        # 命令行
        cmd = 'module purge&&module load vina&&' + cmd
        # 执行
        os.system(cmd)

    def ligand_prep(self, cmd, ligand):
        # 分子名
        lig_name = ligand.split('.')[0]
        # 准备前的分子
        ligand_mol2 = self.ligand_mol2.format(ligand)
        # 准备后的分子
        ligand_pdbqt = self.ligand_pdbqt.format(lig_name)
        # 命令行
        cmd = 'module purge&&module load vina&& timeout 60 ' + cmd.format(ligand_mol2, ligand_pdbqt)
        # 执行
        try:
            os.system(cmd)
        except:
            pass

    def ligand_docking(self, cmd, soft, ligand):
        # 分子名
        lig_name = ligand.split('.')[0]
        # 对接前的分子
        ligand_pdbqt = self.ligand_pdbqt.format(lig_name)
        # 判断分子是否准备成功
        if os.path.exists(ligand_pdbqt):
            # 打分文件
            log_file = self.docking_log.format(lig_name)
            # 对接后分子
            if soft == 'Smina':
                ligand_docked = self.ligand_docked.format(lig_name)
                # 命令行
                cmd = cmd.format(self.protein_pdbqt, ligand_pdbqt, ligand_docked, log_file)
            else:  # qvina  对接后转换格式
                ligand_docked = ligand_pdbqt
                # 命令行
                cmd = cmd.format(self.protein_pdbqt, ligand_pdbqt, ligand_docked,
                                 log_file) + '&& module load openbabel&&obabel {} -O {}'.format(ligand_docked,
                                                                                                self.ligand_docked.format(
                                                                                                    lig_name))
            # 执行
            try:
                os.system(cmd)
            except:
                pass

    # 获取打分的分数
    def get_score(self, ligand):
        # 分子名
        lig_name = ligand.split('.')[0]
        # 打分文件
        log_file = self.docking_log.format(lig_name)
        # 若文件存在
        if os.path.exists(log_file):
            # 读取数据
            with open(log_file, 'r') as f:
                con = f.readlines()
            # 整理数据
            for i in range(len(con)):
                if con[i].startswith('-----'):
                    try:
                        result = [lig_name, con[i + 1].split()[1].strip()]  # 获取数据
                        pd.DataFrame(result).T.to_csv(self.docking_score_csv, index=False, header=False,
                                                      mode='a')  # 写数据到csv中
                    finally:
                        break
