#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-10

import sys
import csv
import oddt
from oddt import fingerprints
import numpy as np


def sparse_to_dense(ifp, size, count_bits=True):
    """Converts sparse fingerprint (indices of 'on' bits) to dense (vector of
    ints).

    Parameters
    ----------
    fp : array-like
        Fingerprint on indices. Can be dupplicated for count vectors.

    size : int
        The size of a final fingerprint.

    count_bits : bool (default=True)
        Should the output fingerprint be a count or boolean vector. If `True`
        the dtype of output is `np.uint8`, otherwise it is bool.


    Returns
    -------
    fp : np.array  (shape=[1, size])
        Dense fingerprint in form of a np.array vector.
    """
    ifp = np.asarray(ifp, dtype=np.uint64)
    if ifp.ndim > 1:
        raise ValueError("Input fingerprint must be a vector (1D)")
    sparsed_fp = np.zeros(size, dtype=np.uint8 if count_bits else bool)
    np.add.at(sparsed_fp, ifp, 1)
    return sparsed_fp


# 获取蛋白分子
protein_file = sys.argv[1]
# 获取小分子
ligand_file = sys.argv[2]
# 获取结果csv
csv_file = sys.argv[3]
# 获取指纹类型
finger_type = sys.argv[4]
# 计算小分子纯名
lig_name = ligand_file.split('/')[-1].split('.')[0]
# 读取蛋白分子
protein = next(oddt.toolkit.readfile('pdb', protein_file))
protein.protein = True
# 读取小分子
ligand = next(oddt.toolkit.readfile('mol2', ligand_file))
# 计算FP
fp = None  # 初始化指纹为None
# 判断指纹类型并计算
if finger_type == 'ifp':  # 根据蛋白变化，指纹长度也会改变
    fp = fingerprints.InteractionFingerprint(ligand, protein, strict=True).tolist()
    fp = sparse_to_dense(fp, size=1024).tolist()  # 转换指纹长度至1024
elif finger_type == 'plec':
    fp = fingerprints.PLEC(ligand, protein, depth_ligand=1, depth_protein=5, distance_cutoff=4.5, size=1024,
                         count_bits=True, sparse=False, ignore_hoh=True).tolist()
elif finger_type == 'splif':  # 根据蛋白变化，指纹长度也会改变 没有定义具体的作用类型，相反的，参考了ECFP，设定了作用半径，有点类似于图
    fp = fingerprints.SPLIF(ligand, protein, depth=1, size=1024, distance_cutoff=4.5)
    fp = sparse_to_dense(fp, size=1024).tolist()  # 转换指纹长度至1024
else:  # silirid 把ifp的指纹按照20种标准氨基酸以及cofactor合并而来，共21类，每类8项相互作用，共计21*8 = 168位
    fp = fingerprints.SimpleInteractionFingerprint(ligand, protein, strict=True).tolist()
# 合并列表
result = [lig_name] + fp
# 把结果写入csv
mycsv = open(csv_file, 'a')
mycsvwriter = csv.writer(mycsv)
mycsvwriter.writerow(result)
mycsv.close()
