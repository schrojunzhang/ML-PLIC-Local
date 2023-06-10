#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-22
# used for Maximum Common Substructure generation

import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from multiprocessing import Pool


def read_mol(ligand_file):
    try:
        mol = Chem.MolFromMolFile(ligand_file)
    except:
        mol = None
        print('read mol failed {}'.format(ligand_file.split('/')[-1].split('.')[0]))
    return mol


def get_fps(mol):
    try:
        # 计算摩根指纹
        fps = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
        # 返回
        return fps
    except:
        print('get ecfp failed')


def cluster_fps(fps, cut_off=0.2):
    # first generate the distance matrix:

    dists = []

    nfps = len(fps)
 
    for i in range(1, nfps):

        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])

        dists.extend([1 - x for x in sims]) 

    # now cluster the data:

    cs = Butina.ClusterData(dists, nfps, cut_off, isDistData=True)

    return cs


def main():
    # 获取复合物名字
    complexes = os.listdir(path_for_complex)
    # mol小分子路径
    ligands = ['{}/{}.mol'.format(path_for_dock, i.split('.')[0]) for i in complexes if
               os.path.exists('{}/{}.mol'.format(path_for_dock, i.split('.')[0]))]
    # 多进程读取数据，计算指纹
    pool = Pool(28)
    mols = pool.map(read_mol, ligands)
    fps = pool.map(get_fps, mols)
    pool.close()
    pool.join()
    # filter None in the fps
    fps = list(filter(None, fps))
    # 计算距离并返回聚类
    clusters = cluster_fps(fps, cut_off=0.85)
    # 获取每类的第一个分子对象
    cluster_first = [mols[cluster[0]] for cluster in clusters]
    # 计算二维坐标，否则生成的图是三维图
    tmp = [AllChem.Compute2DCoords(m) for m in cluster_first]
    # 绘图
    img = Draw.MolsToGridImage(cluster_first, molsPerRow=3, subImgSize=(382, 300),
                           legends=['cluster_{}'.format(i) for i in range(1, len(cluster_first) + 1)])
    # 输出图片
    img.save(cluster_png)



'''
  # # 计算最大子结构
  # res = rdFMCS.FindMCS(mols)
  # # 获取smarts
  # smarts = res.smartsString
  # # 获取分子对象
  # mol = Chem.MolFromSmarts(smarts)
  # # 作图
  # Draw.MolToFile(mol, mcs_png)
  # # 输出最大子结构
  # pd.DataFrame(['smarts', [smarts]]).to_csv(data_csv, index=False, header=False)
'''

if __name__ == '__main__':
    path_for_complex = sys.argv[1]  # 复合物路径
    path_for_dock = sys.argv[2]  # mol文件路径
    # mcs_png = sys.argv[3]  # 最大子结构图
    cluster_png = sys.argv[3]  # 聚类结构图
    # data_csv = sys.argv[4]  # 存放smarts的csv
    main()
