#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-10

import sys
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
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


# 获取小分子
ligand_file = sys.argv[1]
# 获取结果csv
csv_file = sys.argv[2]
# 获取指纹类型
finger_type = sys.argv[3]
# 计算小分子纯名
lig_name = ligand_file.split('/')[-1].split('.')[0]
# 读取小分子
ligand = Chem.MolFromMolFile(ligand_file)
# 计算FP
fp = None  # 初始化指纹为None
# 判断指纹类型并计算
if finger_type == 'ecfp':
    fp = AllChem.GetMorganFingerprintAsBitVect(ligand, 2, nBits=1024)
    fp = fp.ToBitString()
else:
    fp = None
# 合并列表
result = [lig_name] + [i for i in fp]
# 把结果写入csv
mycsv = open(csv_file, 'a')
mycsvwriter = csv.writer(mycsv)
mycsvwriter.writerow(result)
mycsv.close()
