#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-23
import os
import sys
import pandas as pd
from multiprocessing import Pool
from base_scripts.base_class import interaction_base


def main():
    # 实例化对象
    this_interaction = interaction_base.oddt_interaction(job_id)
    # 定义函数
    def each_lig_interaction(lig_names, out_lis=False):
        ligands = ['{}/{}.mol'.format(this_interaction.path_for_dock, lig_name) for lig_name in lig_names if
                   os.path.exists('{}/{}.mol'.format(this_interaction.path_for_dock, lig_name))]
        # 组成列表
        length = len(ligands)
        ligands_hb = zip(ligands, ['hb']*length, [out_lis]*length)
        ligands_clb = zip(ligands, ['clb']*length, [out_lis]*length)
        ligands_qq = zip(ligands, ['qq']*length, [out_lis]*length)
        ligands_lipo = zip(ligands, ['lipo']*length, [out_lis]*length)
        ligands_metal = zip(ligands, ['metal']*length, [out_lis]*length)
        # 多进程计算 计算H键作用
        pool = Pool(28)
        print('cal Hb')
        hb_infos = pool.starmap(this_interaction.interaction2recNum, ligands_hb)
        print('cal clb')
        clb_infos = pool.starmap(this_interaction.interaction2recNum, ligands_clb)
        print('cal qq')
        qq_infos = pool.starmap(this_interaction.interaction2recNum, ligands_qq)
        print('cal lipo')
        lipo_infos = pool.starmap(this_interaction.interaction2recNum, ligands_lipo)
        print('cal metal')
        metal_infos = pool.starmap(this_interaction.interaction2recNum, ligands_metal)
        pool.close()
        pool.join()
        # 合并数据
        interaction_infos = zip(lig_names, hb_infos, clb_infos, qq_infos, lipo_infos, metal_infos)
        # 写入CSV
        pd.DataFrame(interaction_infos, columns=['name', 'hb', 'halogenbond', 'salt_bridges', 'hydrophobic', 'metal']).to_csv(this_interaction.interaction_csv, index=False)
    # 通过描述符文件是否存在判断是对接任务还是pipeline任务
    if os.path.exists(this_interaction.descriptor_csv):  # pipeline
        df = pd.read_csv(this_interaction.descriptor_csv, encoding='utf-8')
        # 获取活性分子
        lig_names = df[df.iloc[:, -1] == 1].iloc[:, 0].values
    else:  # 对接模式
        # 获取小分子
        lig_names = [i.split('.')[0] for i in os.listdir(this_interaction.path_for_complex)]
    # 计算每个分子的相互作用
    each_lig_interaction(lig_names)
    # 计算频次,写入CSV  pipeline模式下，相互作用的统计用训练集中的活性分子的数据
    pool = Pool(28)
    pool.map(this_interaction.get_rec_frequence, this_interaction.interactions.keys())
    pool.close()
    pool.join()
    # 获取被预测为活性的小分子
    lig_names = [i.split('.')[0] for i in os.listdir(this_interaction.path_for_complex)]
    # 计算预测出的化合物的相互作用
    each_lig_interaction(lig_names, out_lis=True)


if __name__ == '__main__':
    job_id = sys.argv[1]  # 获取任务名
    main()
