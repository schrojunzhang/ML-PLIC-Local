#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-09

import sys
from NNScore2module import PDB, binana, command_line_parameters

# 蛋白
protein = sys.argv[1]
# 小分子
ligand = sys.argv[2]
# 日志文件
log_file = sys.argv[3]
# 指令
cmd = "/home/xujun/Soft/SCORE_Function/NNscore/NNScore2module.py -receptor {0} -ligand {1}".format(protein, ligand)
# 初始化空列表
result = []
# 执行
try:
    params_list = cmd.split()
    cmd_params = command_line_parameters(params_list)
    receptor = PDB()
    receptor.LoadPDB_from_file(protein)
    receptor.OrigFileName = protein
    d = binana(ligand, receptor, cmd_params, "", "", "")
    result = d.vina_output + d.ligand_receptor_atom_type_pairs_less_than_two_half.values() + d.ligand_receptor_atom_type_pairs_less_than_four.values() \
                     + d.ligand_atom_types.values() + d.ligand_receptor_atom_type_pairs_electrostatic.values() + d.rotateable_bonds_count.values() \
                     + d.active_site_flexibility.values() + d.hbonds.values() + d.hydrophobics.values() + d.stacking.values() + d.pi_cation.values() \
                     + d.t_shaped.values() + d.salt_bridges.values()
finally:
    # 写TXT文件
    with open(log_file, 'w')as f:
        f.write(str(result))
