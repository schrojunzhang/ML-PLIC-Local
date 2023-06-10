#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-16
# used for cal interaction

import sys
import pandas as pd
import numpy as np
from oddt.toolkits.rdk import readfile
from oddt.interactions import hbond_acceptor_donor, halogenbond_acceptor_halogen, salt_bridges, hydrophobic_contacts, acceptor_metal
from base_scripts.base_class.global_base import job_control

class oddt_interaction(job_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 蛋白质mol对象
        self.protein_mol = readfile(format='pdb', filename=self.protein).__next__()
        # 蛋白质文本对象
        self.protein_data = self.get_protein_data()
        # 金属元素类型
        self.metal_atoms = ['FE', 'MG', 'MN', 'ZN', 'CA', 'NA']
        # 不同的相互作用
        self.interactions = {
            'hb': hbond_acceptor_donor,  # 有cutoff参数，可以设置判定的阈值
            'qq': salt_bridges,
            'clb': halogenbond_acceptor_halogen,
            'lipo': hydrophobic_contacts,
            'metal': acceptor_metal
        }
        # 不同相互作用到列的映射
        self.interactions2col = {
            'hb': 1,  # 有cutoff参数，可以设置判定的阈值
            'clb': 2,
            'qq': 3,
            'lipo': 4,
            'metal': 5
        }

    # 获取配体mol对象
    def get_mol(self, ligand_file):
        return readfile(format='mol', filename=ligand_file).__next__()

    # 获取蛋白质文件字符串信息
    def get_protein_data(self):
        # 打开pdb文件
        with open(self.protein, 'r') as f:
            protein_data = f.readlines()
        # 返回
        return protein_data

    # 获取残基信息
    def get_rec(self, atm_num):
        # 按行读取数据
        for i in range(atm_num, len(self.protein_data)):
            # 获取对应行数据
            if str(atm_num + 1) in self.protein_data[i]:
                # 对行数据按空格分割
                rec_data = self.protein_data[i].split()
                # 残基元素
                rec_element = rec_data[-1]
                # 残基名称
                rec_name = rec_data[3]
                # 残基编号
                rec_id = rec_data[4]
                return str([rec_name, rec_id, rec_element])
        return None


    def has_metal(self, ligandorprotein):
        # 读取数据
        with open(ligandorprotein, 'r') as f:
            content = f.read()
        # 逐行判断
        for metal in self.metal_atoms:
            if metal in content:  # 若蛋白中存在金属元素，返回True
                return True
        # 否则返回False
        return False

    # 计算蛋白配体氢键，并获取蛋白上形成氢键的原子编号
    def interaction2recNum(self, ligand, interaction_type, out_lis):
        # 定义空字典存放残基信息
        rec_infos = {}
        # 获取对象
        mols = [[self.protein_mol, self.get_mol(ligand)], [self.get_mol(ligand), self.protein_mol]]
        # 定义是否只计算一次
        cal_once = False
        # 分别计算 蛋白做受体和供体时的残基情况
        for i in [0, 1]:
            # 若为金属作用
            if interaction_type == 'metal':
                # 若蛋白和小分子都有金属
                if self.has_metal(ligand) and self.has_metal(self.protein):
                    data_array = self.interactions[interaction_type](mols[i][0], mols[i][1])  # 交替次序
                    # 小分子含有金属
                elif self.has_metal(ligand):
                    data_array = self.interactions[interaction_type](mols[0][0], mols[0][1])  # 小分子在后
                    cal_once = True
                    # 蛋白含有金属或者都没有金属
                else:
                    data_array = self.interactions[interaction_type](mols[1][0], mols[1][1])  # 蛋白在后
                    cal_once = True
            # 除金属外的其他作用
            else:
                data_array = self.interactions[interaction_type](mols[i][0], mols[i][1])
            protein_infos = data_array[i]  # 获取供体列表
            # 判断是否含有该相互作用
            if len(protein_infos) == 0:  # 不含相互作用
                pass
            else:  # 含有相互作用
                protein_infos = set([i[0] for i in protein_infos])
                # 获取残基信息
                for atm_num in protein_infos:  # 获取原子的编号
                        atm_num = self.get_rec(atm_num)  # 获得残基信息 [MET, 86, O]
                        rec_infos[atm_num] = rec_infos.get(atm_num, 0) + 1  # 统计次数 dic[None] = {}
                # 判断是否只计算一次
                if cal_once:
                    break
        # 数据样式 {"['MET', '86', 'O']": 1, "['MET', '86', 'N']": 1}
        # 判断是否整理成list的格式
        if out_lis:
            # 空列表暂存
            tmp_lis = []
            # 循环
            for rec_info, rec_fre in rec_infos.items():
                # 把残基信息由字符串变为字典
                rec_info = eval(rec_info)
                # ， 并加入频率
                rec_info.append(rec_fre)
                # 汇总加入新列表
                tmp_lis.append(rec_info)
            # 判断该相互作用残基是否全为空
            if len(tmp_lis) != 0:
                tmp_lis.sort(key=lambda x: x[-1], reverse=True)  # 排序
            # 变量
            rec_infos = tmp_lis
        # 返回数据
        return rec_infos




    def get_rec_frequence(self, interaction_type):
        # 初始化字典存放变量
        frequencies = {}
        # 读数据
        df = pd.read_csv(self.interaction_csv, encoding='utf-8').dropna()
        # 获取列
        interactions = df.iloc[:, self.interactions2col[interaction_type]].values
        # 取出每个分子的相互作用
        for lig_interaction in interactions:
            lig_interaction = eval(lig_interaction)  # 由字符串变回字典
            # 取出字典中的每个残基
            for rec in lig_interaction:
                frequencies[rec] = lig_interaction.get(rec, 0) + frequencies.get(rec, 0)  # 增加次数
        # 由字典变成列表
        # 获取残基信息和频率
        frequency_lis = []
        for rec_info, rec_fre in frequencies.items():
            # 把残基信息由字符串变为字典
            rec_info = eval(rec_info)
            # ， 并加入频率
            rec_info.append(rec_fre)
            # 汇总加入新列表
            frequency_lis.append(rec_info)
        # 判断该相互作用残基是否全为空
        if len(frequency_lis) != 0:
            frequency_lis.sort(key=lambda x: x[-1], reverse=True)  # 排序
            # 写入csv
            pd.DataFrame([interaction_type, frequency_lis]).T.to_csv(self.data_collect_csv, header=False, index=False, mode='a')
            # con = '{}\n{}\n'.format('hb', frequency_lis)  # 定义文本内容
            # with open(self.data_collect_csv, 'a') as f:  # 写入文件
            #     f.write(con)
        else:  # 若全空，返回空列表
            pd.DataFrame([interaction_type, []]).T.to_csv(self.data_collect_csv, index=False, header=False, mode='a')  # 写入CSV
