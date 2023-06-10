#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-08

import sys
from ccdc.protein import Protein
from ccdc.io import MoleculeWriter

protein_file = sys.argv[1]  # 原蛋白文件
pre_protein_file = sys.argv[2]  # 处理后蛋白文件

mol = Protein.from_file('%s' % protein_file)  # 读取文件
mol.remove_all_waters()  # 移除水分子
mol.remove_unknown_atoms()  # 移除未知原子
mol.add_hydrogens()  # 加H
with MoleculeWriter('%s' % pre_protein_file) as protein_writer:
    protein_writer.write(mol)  # 输出文件