#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-08

import os
import sys
import pybel

temp_ligand = sys.argv[1]  # sdf小分子
dst_path = sys.argv[2]  # 文件路径
lig_name = sys.argv[3]  # 小分子纯名字
dst_ligand = '{0}/{1}.mol2'.format(dst_path, lig_name)  # 加上电荷后的分子
add_py = '{0}/{1}_addcharge.py'.format(dst_path, lig_name)  # 定义加电荷脚本
mol = pybel.readfile("sdf", temp_ligand).next()  # 读取分子
net_charge = mol.charge  # 获取电荷
#  定义脚本内容
con = '''from chimera import runCommand
runCommand("open 2 {0}")
runCommand("addcharge nonstd #2 {1} method am1")
runCommand("write format mol2 atomTypes sybyl 2 {2}")
runCommand("close 2")
'''.format(temp_ligand, net_charge, dst_ligand)
# 创建文件
with open(add_py, 'w') as f:
    f.write(con)
# 执行加电荷操作
cmdline = 'cd {0} &&'.format(dst_path)
cmdline += 'module load chimera &&'
cmdline += 'chimera --nogui --nostatus --script {0} &&'.format(add_py)
cmdline += 'rm {0} '.format(add_py)
os.system(cmdline)
