#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-04

import os
import random
import time
import glob
import re
import shutil
import pandas as pd
import numpy as np
from Crypto.Cipher import AES
from Crypto.Util.Padding import pad
from Crypto.Util.Padding import unpad
from django.core.mail import EmailMultiAlternatives
from multiprocessing import Pool


class job_generation(object):
    def __init__(self):
        self.file_path = '/home/xujun/Project_1/mlplic/submitted_files'  # 定义文件存放路径
        self.path_for_download = '/home/xujun/Project_1/mlplic/files_for_download'  # example文件下载路径
        # 定义密钥
        self.__pwd = b'xxxx'
        self.__iv = b'xxxx'
        # 定义密钥文件
        self.key_file = self.file_path + '/{}/config/id.rsa'
        # 访问次数文件
        self.total_viewed = '{}/viewed_time.txt'.format(self.path_for_download)

    # 生成唯一ID
    def get_id(self):
        # 随机生成4位数
        rand = random.randint(1000, 9999)
        # 生成当前时间戳,取后四位
        t_string = str(time.time()).replace('.', '1')[-4:]
        # 初始化id
        job_id = 'job_{}{}'.format(rand, t_string)
        # 构建递归函数
        self.path_job = '{}/{}'.format(self.file_path, job_id)
        # 检查任务名是否存在
        if os.path.exists(self.path_job):
            self.get_id()
        else:
            return job_id

    # 加密解密文件
    def make_secure(self, content, encryption=True):
        # 实例化对象
        cryptor = AES.new(self.__pwd, AES.MODE_CBC, self.__iv)
        # 判断是加密还是解密
        # 加密  content 是 字符串
        if encryption:
            content = cryptor.encrypt(pad(content.encode('utf-8'), AES.block_size))
            return content
        # 解密  content 是 bytes
        else:
            # 判断数字是否可以为16位
            if len(content) % 16 == 0:
                content = unpad(cryptor.decrypt(content), AES.block_size)
                return content.decode('utf-8')
            else:  # 若不是16位，返回None
                return None

    # 获取访问量
    def get_total_viewed(self):
        # 打开访问次数文件
        with open(self.total_viewed, 'r') as f:
            content = f.read()
        # 获取访问次数
        total_viewed = int(content.strip())
        return total_viewed


class job_control(job_generation):

    def __init__(self, job_id):
        super().__init__()
        # 任务类型
        self.job_type = '#descriptors'
        # 全局文件路径
        self.sfpy_path = '/home/xujun/Project_1/mlplic/base_scripts/base_descriptors'  # 定义计算打分函数脚本所在路径
        self.model_path = '/home/xujun/Project_1/mlplic/base_scripts/base_models'  # 定义模型脚本所在路径
        self.docking_path = '/home/xujun/Project_1/mlplic/base_scripts/base_docking'  # 定义对接脚本所在路径
        self.screening_path = '/home/xujun/Project_1/mlplic/base_scripts/base_screening'  # 定义筛选脚本所在路径
        self.help_path = r'/home/xujun/Project_1/mlplic/base_scripts/help_scripts'  # 辅助脚本路径
        self.ifp_module = 'source activate my-rdkit-env\npython'  # 相互作用指纹执行前环境设置
        # 定义html的<input name = 'sf'>与打分函数全名的映射
        self.scoring_function = {
            'aff': 'affiscore',
            'affi': 'affinitydG',
            'alp': 'alphaHB',
            'ase': 'ase',
            'asp': 'asp',
            'aut': 'autodock',
            'gau': 'chemgauss',
            'che': 'chemplp',
            'chem': 'chemscore',
            'con': 'contact',
            'cont': 'continuous',
            'cys': 'cyscore',
            'dsx': 'dsx',
            'gal': 'galaxy',
            'gbv': 'GBVIWSAdG',
            'sp': 'sp',
            'xp': 'xp',
            'gol': 'goldscore',
            'gri': 'grid',
            'haw': 'hawkins',
            'ldg': 'LondondG',
            'nn': 'nnscore',
            'plp': 'plp',
            'plp9': 'plp95',
            'rdo': 'rdock',
            'rdoc': 'rdock_sol',
            'cre': 'rfscore_credo',
            'ele': 'rfscore_element',
            'syb': 'rfscore_sybyl',
            'sas': 'sasa',
            'smi': 'smina',
            'smo': 'smog2016',
            'sur': 'surflex',
            'vin': 'vina',
            'xsc': 'xscore',
            'ifp': 'ifp',
            'ple': 'plec',
            'sil': 'silirid',
            'spl': 'splif',
            'ecf': 'ecfp',
            'pad': 'pubchem'}
        # 定义html的<input name = ''>与值的映射  口袋坐标
        self.binding_site = {'site_xyz': self.get_pocket_from_xyz, 'site_res': self.get_pocket_from_res}
        # 上传文件格式映射
        self.file_extension = {'protein_file': 'pdb',
                               'ligand_file': 'mol2',
                               'crystal_file': 'mol2',
                               'decoys_file': 'mol2',
                               'test_file': 'mol2',
                               'descriptor_file': 'csv',
                               'test_descriptor_file': 'csv',
                               'label_file': 'csv'}
        # 定义具体任务各级路径
        self.path_job = '{}/{}'.format(self.file_path, job_id)
        self.path_for_dock = '{}/undocked'.format(self.path_job)  # 用于存放未对接的原始小分子
        self.path_for_lig = '{}/files'.format(self.path_job)  # 蛋白小分子存放路径
        self.path_for_log = '{}/log'.format(self.path_job)  # 定义日志路径
        self.path_for_csv = '{}/csvs'.format(self.path_job)  # 定义描述符存放路径
        self.path_for_config = '{}/config'.format(self.path_job)  # 定义配置存放路径
        self.path_for_complex = '{}/complex'.format(self.path_job)  # 定义存放用于查看PL相互作用的蛋白配体复合物
        self.path_for_report = '{}/report'.format(self.path_job)  # 报告文件存放路径
        # 定义具体文件
        # 分子文件
        self.protein = '{}/protein_file.pdb'.format(self.path_for_lig)  # 上传的蛋白分子文件
        self.ligand = '{}/ligand_file.mol2'.format(self.path_for_lig)  # 上传的配体文件
        self.decoys = '{}/decoys_file.mol2'.format(self.path_for_lig)  # 上传的无活性分子文件
        self.test_ligand = '{}/test_file.mol2'.format(self.path_for_lig)  # 上传的待测分子
        self.ligands = [i.split('/')[-1] for i in glob.glob('{}/*.mol2'.format(self.path_for_lig))]  # 分割后的小分子名字集
        self.crystal_ligand = '{}/crystal_file.mol2'.format(self.path_for_lig)  # 上传或生成的晶体结构
        self.config_file = '{}/config'.format(self.path_for_config)  # 定义任务配置文件
        self.complex_file = self.path_for_complex + '/{}.pdb'  # 复合物分子
        self.pocket_file = '{}/pocket.pdb'.format(self.path_for_report)  # 口袋分子
        # 脚本文件
        self.split_py = '{}/ligand_split.py'.format(self.help_path)  # 分割分子脚本路径
        self.merge_py = '{}/merge_csv.py'.format(self.help_path)  # 合并CSV脚本路径
        self.docking_py = '{}/docking.py'.format(self.docking_path)  # 对接脚本路径
        self.screening_py = '{}/screening.py'.format(self.screening_path)  # 虚筛脚本路径
        self.result_deal_py = '{}/deal_result.py'.format(self.help_path)  # 结果处理脚本路径
        # csv文件
        self.descriptor_csv = '{}/descriptor_file.csv'.format(self.path_for_csv)  # 汇总的csv文件
        self.test_descriptors_csv = '{}/test_descriptor_file.csv'.format(self.path_for_csv)  # 汇总的测试集特征csv
        self.docking_score_csv = '{}/score.csv'.format(self.path_for_report)  # 对接打分值
        self.data_collect_csv = '{}/data.csv'.format(self.path_for_report)  # smarts及残基频率分析数据的存放CSV
        self.bias_var_csv = '{}/bias_var.csv'.format(self.path_for_report)  # 方差偏差分析数据存放的CSV
        self.interaction_csv = '{}/interaction.csv'.format(self.path_for_report)  # 各种相互作用的存放CSV
        self.model_metric_csv = '{}/model_metric.csv'.format(self.path_for_report)  # 模型评价指标的存放CSV
        self.feature_importance_csv = '{}/feature_importance.csv'.format(self.path_for_report)  # 特征重要性的存放CSV
        self.final_parameter_csv = '{}/final_parameter.csv'.format(self.path_for_report)  # 特征重要性的存放CSV
        self.active_ligand_csv = '{}/active_ligand.csv'.format(self.path_for_report)  # 被预测为抑制剂的小分子名字的存放CSV
        # 图片
        self.bias_var_png = '{}/bias_var.png'.format(self.path_for_report)  # 方差偏差分析图
        self.learning_curve_png = '{}/learning_curve.png'.format(self.path_for_report)  # 学习曲线分析图
        self.hyperParameter_scatter_png = '{}/hyperParameter_scatter.png'.format(self.path_for_report)  # 超参数寻优图
        self.roc_png = '{}/roc.png'.format(self.path_for_report)  # ROC面积寻优图
        self.shap_png = '{}/shap.png'.format(self.path_for_report)  # shap模型解释图
        # 结果文件
        self.result_file = '{}/result.zip'.format(self.path_job)  # 结果压缩文件
        # pbs文件 用于判断是否发送过邮件
        self.pbs_file = '{}/{}_submit.pbs'.format(self.path_for_log, job_id)
        # 邮件文字内容
        # 任务提交成功
        self.email_submit_content = '{} has been successfully submitted.The key file used for getting access to the resu' \
                                    'lt is in the attachment.Once the job finished, you will be informed by email.'.format(
            job_id)
        # 任务完成
        self.job_finished_content = '{} has finished.You can get the result directly through the link below or go to the ' \
                                    'queue interface to get the result.http://cadd.zju.edu.cn/mlplic/results/{}'.format(
            job_id, self.make_secure(job_id, encryption=True))

    def handle_uploaded_file(self, request, input_name, dest_path):  # 上传文件
        if input_name in request.FILES:  # 判断是否含有该文件数据
            src_file = request.FILES[input_name]  # 待上传文件
            # 判断文件是否为config
            if input_name == 'config_file':
                dst_file = self.config_file  # 若为config，则重定位文件路径
            elif input_name == 'descriptor_file':  # 若为descriptor
                dst_file = self.descriptor_csv
            elif input_name == 'test_descriptor_file':
                dst_file = self.test_descriptors_csv
            else:  # 若不是config，
                file_extension = self.file_extension  # 获取文件后缀
                dst_file = '{}/{}.{}'.format(dest_path, input_name, file_extension[input_name])  # 上传到服务器的路径
            with open(dst_file, 'wb+') as destination:
                for chunk in src_file.chunks():  # 分片段上传
                    destination.write(chunk)

    def check_name(self, lig_name, i, dst_path):  # 检查名字的独一性
        # 正则表达式判断获取的名字是否有误
        name_modern = re.compile(r'^[A-Za-z0-9]+[0-9a-zA-Z\_]*[A-Za-z0-9]*$')  # 正则表达式 由英文字母和数字组成
        if not name_modern.match(lig_name):
            lig_name = i  # 若有误，则用分子序号当作名字
        lig_file = '{}/{}.mol2'.format(dst_path, lig_name)  # 定义输出分子路径

        # 递归检查重名在文件
        def same_file(lig_file, n=0):
            if os.path.exists(lig_file):
                n += 1
                lig_file = '{}/{}_{}.mol2'.format(dst_path, lig_name, n)
                return same_file(lig_file, n)
            else:
                return lig_file

        lig_file = same_file(lig_file)  # 检查是否重名
        return lig_file

    def split_file(self, src_ligand, dst_path, return_name=False):  # 分割小分子(mol2)
        if not return_name:  # 不加标签
            # 读取数据，存到con中
            with open(src_ligand, 'r') as f:
                con = f.read()
            # 根据@<TRIPOS>MOLECULE分割字符串
            con = con.split('@<TRIPOS>MOLECULE')
            for i in range(1, len(con)):
                # lig_name = con[i].split('\n')[1].strip()  # 获取小分子名字
                lig_name = con[i].splitlines()[1].strip()  # 获取小分子名字
                # 检查名字的独一性
                lig_file = self.check_name(lig_name, i, dst_path)
                # 创建文件
                with open(lig_file, 'w') as f:
                    f.write('@<TRIPOS>MOLECULE' + con[i])
            # 分裂结束后删除原小分子
            os.remove(src_ligand)
            # 判断是否存在晶体结构
            crystal_file = self.crystal_ligand
            if not os.path.exists(crystal_file):
                shutil.copyfile(lig_file, crystal_file)  # 把最后一个分子当作晶体结构
        else:  # 加标签
            lig_names = []  # 初始化列表存放名字
            # 读取数据，存到con中
            with open(src_ligand, 'r') as f:
                con = f.read()
            # 根据@<TRIPOS>MOLECULE分割字符串
            con = con.split('@<TRIPOS>MOLECULE')
            for i in range(1, len(con)):
                lig_name = con[i].split('\n')[1].strip()  # 获取小分子名字
                # 检查名字的独一性
                lig_file = self.check_name(lig_name, i, dst_path)
                # 把小分子纯名字加入列表
                lig_names.append(lig_file.split('/')[-1].split('.')[0])
                # 创建文件
                with open(lig_file, 'w') as f:
                    f.write('@<TRIPOS>MOLECULE' + con[i])
            # 分裂结束后删除原小分子
            os.remove(src_ligand)
            # 若分子为活性分子
            if src_ligand == self.ligand:
                # 判断是否存在晶体结构
                crystal_file = self.crystal_ligand
                if not os.path.exists(crystal_file):
                    shutil.copyfile(lig_file, crystal_file)  # 把活性分子的最后一个分子当作晶体结构
            # 返回名字列表
            return lig_names

    # 用于descriptors 和 pipline 不对接的情况下
    def split_addLabel(self):  # 切割文件并加上标签
        # 分别对decoys和ligand进行分割和加标签
        if not os.path.exists(self.decoys):  # 若decoys及test文件不存在，则只进行分割文件
            self.split_file(src_ligand=self.ligand, dst_path=self.path_for_lig, return_name=False)
        else:
            # 若上传test分子， 生成test_descriptor.csv
            if os.path.exists(self.test_ligand):
                # 获取名字
                test_names = self.split_file(self.test_ligand, dst_path=self.path_for_lig, return_name=True)
                # 写入test_descriptor_csv
                pd.DataFrame(test_names, columns=['name']).to_csv(self.test_descriptors_csv, index=False)
            # 不管是否上传了test分子，获取decoys的名字并标记
            active_names = ['crystal_file'] + self.split_file(self.ligand, dst_path=self.path_for_lig,
                                                              return_name=True)  # 加上单独的crystalfile
            decoys_names = self.split_file(self.decoys, dst_path=self.path_for_lig, return_name=True)
            active_names.extend(decoys_names)  # 合并名字
            lig_names = active_names  # 把全名赋给变量lig_names
            # 赋予标签
            labels = [1 for i in range(len(active_names) - len(decoys_names))] + [0 for i in range(len(decoys_names))]
            # 写入descriptors.csv
            df = pd.DataFrame(lig_names, columns=['name'])  # 写入名字
            df['class'] = pd.DataFrame(labels)  # 写入标签
            df.to_csv(self.descriptor_csv, index=False)  # 创建csv

    # 从config里读取表单
    def read_config(self, job_type):
        # 初始化表单
        job_form = None
        # 读取配置文件
        with open(self.config_file, 'rb') as f:
            content = f.read()  # 读取表单
        # 解密
        content = self.make_secure(content, encryption=False)
        # 分割成行
        content = content.split('\n')
        # 获取对应表单
        for line in content:
            if line.startswith(job_type):  # 寻找标题
                job_form = eval(content[content.index(line) + 1])  # 找到后获取标题下一行内容即为表单
                break
        return job_form

    # 创建config文件
    def write_config(self, content, *, xyz=''):
        # 判断是否输入xyz
        if xyz:  # 若输入xyz,则向原文件中加入口袋信息
            content = xyz
        # 判断当前是否存在config
        if os.path.exists(self.config_file):  # 若存在
            with open(self.config_file, 'rb') as f:
                content = self.make_secure(f.read(), encryption=False) + content  # 把原来文本内容解密,并把新内容加在其后
        # 加密
        content = self.make_secure(content, encryption=True)
        # 写文件
        with open(self.config_file, 'wb') as f:
            f.write(content)

    # 用pandas合并csv
    def merge_csv(self, dst_file):
        # 获得CSV文件
        csvs = os.listdir(self.path_for_csv)
        # 判断训练集描述符文件是否存在
        # 单纯计算描述符的情况
        if not os.path.exists(self.descriptor_csv) and not os.path.exists(self.test_descriptors_csv):
            df_total = pd.read_csv('{}/{}'.format(self.path_for_csv, csvs[0]))  # 读取第一个文件
            csvs.remove(csvs[0])  # 全体文件列表去除第一份文件
        # 上传了decoys
        elif not os.path.exists(self.test_descriptors_csv):
            # 通过互斥csv的名字
            uniqe_name = self.descriptor_csv.split('/')[-1]
            # 移除互斥csv
            csvs.remove(uniqe_name)
            # 读取第一份数据
            df_total = pd.read_csv('{}'.format(dst_file))
        # 分割训练集和测试集 （即用户上传了decoys和testfile）
        else:
            # 通过互斥csv的名字
            uniqe_name = [self.descriptor_csv.split('/')[-1], self.test_descriptors_csv.split('/')[-1]]
            # 移除互斥csv
            for csv in uniqe_name:
                if csv in csvs:
                    csvs.remove(csv)
            # 读取第一份数据
            df_total = pd.read_csv('{}'.format(dst_file))
        # 统一name的数据类型
        df_total.name = df_total.name.astype(str)
        # 合并剩下所有的csv
        for i in range(0, len(csvs)):
            # 文件
            csv_file = '{}/{}'.format(self.path_for_csv, csvs[i])
            # df
            df = pd.read_csv(csv_file, encoding='utf-8')
            # 统一name数据类型
            df.name = df.name.astype(str)
            # 合并数据
            df_total = pd.merge(left=df_total, right=df, how='left', on='name')
        # 获取列名
        columns = df_total.columns.tolist()
        # 判断 标签是否在df中
        if 'class' in columns:
            columns.remove('class')  # 原列名中移除class
            columns += ['class']  # class加到最后一列
            df_total = df_total[columns]  # 转换df
        # 输出文件
        df_total.to_csv(dst_file, index=False)

    # 残基获取口袋
    def get_pocket_from_res(self, xyz):
        binding_ress = xyz
        # 按行读取蛋白数据
        with open(self.protein, 'r') as f:
            pdb_data = f.readlines()
        # 获取残基号
        # 'HIS_14, THR_83, PHE_150'
        binding_ress = binding_ress.split(',')
        # ['HIS_14', ' THR_83', ' PHE_150']
        binding_ress = [binding_res.strip() for binding_res in binding_ress]
        # ['HIS_14', 'THR_83', 'PHE_150']
        # 获取残基号
        binding_ress = [float(binding_res.split('_')[-1]) for binding_res in binding_ress]
        # 排序
        binding_ress.sort()
        # 获取残基号
        xyz = []  # 坐标列表
        max_res = binding_ress[-1]  # 最大残基号
        min_res = binding_ress[0]  # 最小残基号
        for res_num in binding_ress:
            # 取出数据
            for line in pdb_data:
                try:
                    # pdb中的残基号
                    pdb_res_num = float(line[23:26].strip())
                    # 若残基号正确
                    if pdb_res_num == res_num:
                        x = float(line[31:38].strip())
                        y = float(line[39:46].strip())
                        z = float(line[47:54].strip())
                        xyz.append([x, y, z])
                except:
                    pass
        # 转成np.array
        xyz = np.array(xyz)
        # 计算均值
        xyz = xyz.mean(axis=0).tolist()
        # 获取口袋信息
        atom_num = 0  # 原子编号
        res_data = '''MODEL        1\n'''  # 残基信息
        for i in range(len(pdb_data)):
            # 每一行的数据
            line = pdb_data[i]
            try:
                # pdb中的残基号
                pdb_res_num = float(line[23:26].strip())
                # 判断pdb残基号是否在大小残基号内
                if pdb_res_num >= min_res and pdb_res_num <= max_res:
                    # 原子编号+1
                    atom_num += 1
                    # 写信息
                    res_data += 'ATOM   {:>4}  {}'.format(atom_num, line[13:])
            except:
                pass
        # 加上结束符
        res_data += '''ENDMDL\nEND'''
        # 写口袋文件
        with open(self.pocket_file, 'w') as f:
            f.write(res_data)
        return xyz

    # 从crystal_ligand中获取口袋坐标
    def get_pocket_from_crystal(self):
        # 配体
        crystal_file = self.crystal_ligand
        # 根据共晶配体确定坐标
        x = os.popen(
            "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $3}' | awk '{x+=$1} END {print x/(NR-2)}'" % crystal_file).read()
        y = os.popen(
            "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $4}' | awk '{y+=$1} END {print y/(NR-2)}'" % crystal_file).read()
        z = os.popen(
            "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $5}' | awk '{z+=$1} END {print z/(NR-2)}'" % crystal_file).read()
        xyz = [float(x.strip()), float(y.strip()), float(z.strip())]
        return xyz

    # 从输入的坐标直接获取口袋坐标
    def get_pocket_from_xyz(self, xyz):
        xyz = xyz.split(',')
        x = float(xyz[0].strip())
        y = float(xyz[1].strip())
        z = float(xyz[2].strip())
        return [x, y, z]

    # 记录口袋坐标
    def get_pocket(self):
        # job_form
        print(self.job_type)
        job_form = self.read_config(job_type=self.job_type)
        # xyz
        xyz = None
        # 判断是否从残基或者坐标获取口袋
        for m in self.binding_site:
            if job_form.get(m):
                xyz = self.binding_site[m](job_form[m])
        # 若以上两种方法都不是
        if not xyz:
            xyz = self.get_pocket_from_crystal()
        # 写入表单
        self.write_config(job_form, xyz='#xyz\n{}\n'.format(xyz))
        return xyz

    # 读取口袋坐标
    def read_pocket(self):
        xyz = self.read_config(job_type='#xyz')
        return xyz[0], xyz[1], xyz[2]

    # 在sh中加入各种任务都需要执行的新脚本
    def add_job(self, con, id_name, job_type, request):
        # 加入新内容
        con += '\npython {} {} {}'.format(self.result_deal_py, id_name, job_type)
        # 返回
        return con

    # 提交任务
    def qsub(self, path_log, id_name, sh_file):
        # 提交任务
        cmd = '''cd {}&&pbs_submit.sh {} cu02 28 bash {}'''.format(path_log, id_name, sh_file)
        os.system(cmd)

    # 计算描述符提交任务
    def submit_job(self, request, id_name, path_log):
        # 定义sh文件路径
        sh_file = '{}/descriptors.sh'.format(path_log)
        # 定义文件开头, split ligand
        con = '''#!/bin/bash\nexport PYTHONPATH=$PYTHONPATH:/home/xujun/Project_1/mlplic\nmodule load anaconda3/5.1.0\npython {} {} {}'''.format(
            self.split_py, id_name, self.job_type.strip('#'))
        # 把表单记录到文件中传递参数
        job_form = request.POST.dict()  # POST转成字典
        del job_form['csrfmiddlewaretoken']  # 删除csrf键
        form_con = '''{}\n{}\n'''.format(self.job_type, job_form)  # 组装内容
        self.write_config(form_con)  # 写入文件
        # 计算口袋
        self.get_pocket()
        # 循环判断被选择的描述符
        for sf in self.scoring_function:
            if sf in request.POST:  # 若被选中
                # 创建执行命令变量
                new_con = '''\npython {}/{}.py {}'''.format(self.sfpy_path, self.scoring_function[sf], id_name)
                # 像sh中加入执行命令
                con += new_con
                # 创建对应文件夹
                os.mkdir('{}/{}/{}'.format(self.file_path, id_name, self.scoring_function[sf]))
        # 合并文件
        con += '\npython {0} {1}\npython {0} {1} {2}'.format(self.merge_py, id_name, self.test_descriptors_csv)
        # 加入新的共同任务:压缩结果文件
        con = self.add_job(con, id_name, self.job_type.strip('#'), request)
        # 创建sh文件
        with open(sh_file, 'w') as f:
            f.write(con)
        # 提交任务
        self.qsub(path_log, id_name, sh_file)

    # 发送email
    def send_email(self, email_ads, job_id, content, ath_file):
        # 导入django配置
        import os
        os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ML-PLIC-Local.settings")
        # 发送地址
        from_email = 'xxxxx@qq.com'
        # 接受邮箱地址
        email_ads = [email_ads]
        # 主题
        subject = 'ML-PLIC_{}'.format(job_id)
        # 头部
        headers = {'Message-ID': 'foo'}
        # 实例化对象
        msg = EmailMultiAlternatives(subject, content, from_email, email_ads, headers)
        # 增加附件
        if ath_file:  # 若文件不为None
            msg.attach_file(ath_file)
        # 防止对方接收不了html文本
        msg.attach_alternative(content=content, mimetype="text/html")
        # 发送
        msg.send()

    # 利用openbabel 转换格式
    def openbabel_transform(self, src_file, dst_file):
        # 指令
        cmd = 'module load openbabel&&obabel {} -O {}'.format(src_file, dst_file)
        # 执行
        os.system(cmd)

    # 读取含有某个关键字的最后一行，用于生成复合物
    def get_final_index(self, data, key_word='ATOM'):
        for final_index in range(len(data)):
            # 若此行不包含关键字，上一行包含，则输出该行索引
            if key_word not in data[final_index]:
                if key_word in data[final_index - 1]:
                    return final_index
        # 若没有找到符合条件的索引，输出所有行数
        return len(data)

    # 生成复合物
    def generate_complex(self, active_ligand):
        # 定义分子
        active_file = '{}/{}.mol2'.format(self.path_for_lig, active_ligand)  # 活性mol2分子
        ligand_file = '{}/{}.pdb'.format(self.path_for_complex, active_ligand)  # 活性pdb分子
        complex_file = self.complex_file.format(active_ligand)  # 复合物分子
        # 转换配体到pdb
        self.openbabel_transform(src_file=active_file, dst_file=ligand_file)
        # 打开蛋白文件
        with open(self.protein, 'r') as f:
            protein = f.readlines()
        # 获取INDEX
        final_index = self.get_final_index(protein)
        # 读取对应数据
        protein = protein[:final_index]
        # 打开小分子文件
        with open(ligand_file, 'r') as f:
            ligand = f.readlines()
        # 获取索引
        final_index = self.get_final_index(ligand)
        # 读取对应数据
        ligand = ligand[: final_index]
        # 替换'ATOM  '为HETATM 要多加两个空格，否则替换后 原子序号会后移两位导致PDB格式不标准，3Dmoljs识别不了
        ligand = [i.replace('ATOM  ', 'HETATM') if i.startswith('ATOM') else i for i in ligand]
        # 合并字符串
        complex_content = ''.join(protein) + ''.join(ligand)
        # 写复合物
        with open(complex_file, 'w') as f:
            f.write(complex_content)

    # 移动分子并组成复合物
    def active2complex(self, active_ligands):
        # 多线程
        pool = Pool(28)
        pool.map(self.generate_complex, active_ligands)
        pool.close()
        pool.join()


class docking_control(job_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 任务类型
        self.job_type = '#docking'
        # 对接文件夹下的文件
        # 准备前
        self.protein = '{}/protein_file.pdb'.format(self.path_for_dock)  # 上传的蛋白分子文件
        self.ligand = '{}/ligand_file.mol2'.format(self.path_for_dock)  # 上传的配体文件
        self.decoys = '{}/decoys_file.mol2'.format(self.path_for_dock)  # 上传的无活性分子文件
        self.test_ligand = '{}/test_file.mol2'.format(self.path_for_dock)  # 上传的待测分子
        self.ligands = [i.split('/')[-1] for i in glob.glob('{}/*.mol2'.format(self.path_for_dock))]  # 分割后的小分子名字集
        self.crystal_ligand = '{}/crystal_file.mol2'.format(self.path_for_dock)  # 上传或生成的晶体结构
        # 准备后
        self.protein_pdbqt = '{}/protein_file.pdbqt'.format(self.path_for_dock)  # 蛋白准备后的文件
        self.ligand_mol2 = self.path_for_dock + '/{}'  # 单个准备前分子
        self.ligand_pdbqt = self.path_for_dock + '/{}.pdbqt'  # 单个准备后分子
        self.crystal_ligand_pdbqt = '{}/crystal_file.pdbqt'.format(self.path_for_dock)  # 准备后晶体结构
        self.ligands_pdbqt = [i.replace('.mol2', '.pdbqt') for i in self.ligands if
                              os.path.exists('{}/{}.pdbqt'.format(self.path_for_dock, i.split('.')[0]))]  # 准备成功的分子集合
        # 对接后
        self.ligand_docked = self.path_for_lig + '/{}.mol2'  # 单个对接后分子
        self.protein_docked = '{}/protein_file.pdb'.format(self.path_for_lig)  # 上传的蛋白分子文件

    def submit_job(self, request, id_name, path_log):
        # 定义sh文件路径
        sh_file = '{}/docking.sh'.format(path_log)
        # 定义文件开头
        con = '''#!/bin/bash\nexport PYTHONPATH=$PYTHONPATH:/home/xujun/Project_1/mlplic\nmodule load anaconda3/5.1.0'''
        # 创建执行命令变量 分割分子和对接
        new_con = '''\npython {0} {1} {2}\npython {3} {1}'''.format(self.split_py, id_name, self.job_type.strip('#'), self.docking_py)
        # 把表单记录到文件中传递参数
        job_form = request.POST.dict()  # POST转成字典
        del job_form['csrfmiddlewaretoken']  # 删除csrf键
        form_con = '''{}\n{}\n'''.format(self.job_type, job_form)  # 组装内容
        self.write_config(form_con)  # 写入文件
        # 计算口袋
        self.get_pocket()
        # 像sh中加入执行命令
        con += new_con
        # 加入新的共同任务:压缩结果文件
        con = self.add_job(con, id_name, self.job_type.strip('#'), request)
        # 创建sh文件
        with open(sh_file, 'w') as f:
            f.write(con)
        # 提交任务,并传入表单
        self.qsub(path_log, id_name, sh_file)


class model_control(job_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 任务类型
        self.job_type = '#modelling'
        # 从表单获取模型映射
        self.algorithm = {
            'eXtreme Gradient Boosting': 'xgb',
            'Support Vector Machine': 'svm',
            'Random Forest': 'rf'
        }

    def submit_job(self, request, id_name, path_log):
        # 读取表单
        job_form = request.POST.dict()
        # 获取算法
        algo = job_form.get('algorithm')
        # 定义sh文件路径
        sh_file = '{}/modelling.sh'.format(path_log)
        # 定义文件开头
        con = '''#!/bin/bash\nexport PYTHONPATH=$PYTHONPATH:/home/xujun/Project_1/mlplic\nmodule load anaconda3/5.1.0'''
        # 创建执行命令变量
        new_con = '''\npython {}/{}.py {}'''.format(self.model_path, self.algorithm[algo], id_name)
        # 像sh中加入执行命令
        con += new_con
        # 加入新的共同任务:压缩结果文件
        con = self.add_job(con, id_name, self.job_type.strip('#'), request)
        # 创建sh文件
        with open(sh_file, 'w') as f:
            f.write(con)
        # 提交任务,并传入表单
        self.qsub(path_log, id_name, sh_file)


class screening_control(job_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 任务类型
        self.job_type = '#screening'
        # 优化参数前缀
        self.para_prefix = '#best_parameter'
        # 前面几种任务类型
        self.dock_job_type = '#docking'
        self.descriptor_job_type = '#descriptors'
        self.model_job_type = '#modelling'

    def submit_job(self, request, id_name, path_log):
        # 定义sh文件路径
        sh_file = '{}/screening.sh'.format(path_log)
        # 定义文件开头
        con = '''#!/bin/bash\nexport PYTHONPATH=$PYTHONPATH:/home/xujun/Project_1/mlplic\nmodule load anaconda3/5.1.0'''
        # 加入建模预测脚本
        con += '''\npython {} {}'''.format(self.screening_py, id_name)
        # 加入新的共同任务:压缩结果文件
        con = self.add_job(con, id_name, self.job_type.strip('#'), request)
        # 创建sh文件
        with open(sh_file, 'w') as f:
            f.write(con)
        # 提交任务,并传入表单
        self.qsub(path_log, id_name, sh_file)


class pipline_control(job_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        self.job_id = job_id
        self.job_type = '#pipeline'
        self.docking_form = ['protein_repair', 'protein_delchain', 'protein_cleanup', 'ligand_repair', 'ligand_charge',
                             'ligand_cleanup', 'space_size', 'exhaustiveness', 'num_modes', 'my_email', 'whether_dock']
        self.modelling_form = ['algorithm', 'train_size', 'scaler', 'feature_selector', 'feature_selector_para',
                               'max_eval',
                               'n_estimators', 'learning_rate', 'subsample', 'max_depth', 'gamma', 'min_child_weight',
                               'colsample_bytree', 'colsample_bylevel', 'colsample_bynode', 'reg_alpha', 'reg_lambda',
                               'min_sample_split', 'min_samples_leaf', 'max_features', 'min_impurity_decrease', 'C',
                               'gamma', 'n_estimators_method', 'learning_rate_method', 'subsample_method',
                               'max_depth_method', 'gamma_method', 'min_child_weight_method',
                               'colsample_bytree_method', 'colsample_bylevel_method', 'colsample_bynode_method',
                               'reg_alpha_method', 'reg_lambda_method',
                               'min_sample_split_method', 'min_samples_leaf_method', 'max_features_method',
                               'min_impurity_decrease_method', 'C_method',
                               'gamma_method']

    def whether_dock(self):
        # 获取对接字典
        dock_dic = self.read_config(job_type='#docking')
        # 判断是否对接
        whether_dock = dock_dic.get('whether_dock')
        # 不对接
        if whether_dock == 'no':
            return False
        else:  # 对接
            return True

    def whether_modelling(self):
        # 获取对接字典
        modelling_dic = self.read_config(job_type='#modelling')
        # 判断是否对接
        whether_modelling = modelling_dic.get('algorithm')
        # 不对接
        if whether_modelling == 'No Modelling':
            return False
        else:  # 对接
            return True

    def split_files(self, docked='False'):
        # 判断接下来是否对接
        if docked == 'False':  # 进行对接， 小分子在undocked 文件夹进行分割
            this_dock = docking_control(self.job_id)
            # 分别获取名字列表
            active_names = ['crystal_file'] + self.split_file(this_dock.ligand, dst_path=self.path_for_dock,
                                                              return_name=True)  # 加上单独的crystalfile
            decoys_names = self.split_file(this_dock.decoys, dst_path=self.path_for_dock, return_name=True)
            test_names = self.split_file(this_dock.test_ligand, dst_path=self.path_for_dock, return_name=True)
            active_names.extend(decoys_names)  # 合并活性和诱饵分子名字
            lig_names = active_names  # 把全名赋给变量lig_names
            # 赋予标签
            labels = [1 for i in range(len(active_names) - len(decoys_names))] + [0 for i in
                                                                                  range(len(decoys_names))]
            # 写入descriptor_csv
            df = pd.DataFrame(lig_names, columns=['name'])  # 写入名字
            df['class'] = pd.DataFrame(labels)  # 写入标签
            df.to_csv(self.descriptor_csv, index=False)  # 创建csv
            # 写入test_descriptor_csv
            pd.DataFrame(test_names, columns=['name']).to_csv(self.test_descriptors_csv, index=False)
            # 读取写入口袋
            this_dock.get_pocket()
        else:  # 不进行对接直接分割分子  小分子在files文件夹进行分割
            # 分割训练分子并写入CSV
            self.split_addLabel()
            # 写入口袋
            self.get_pocket()

    def write_config(self, content, *, xyz=''):
        # 初始化表单
        form = '''#docking\n{0}\n#descriptors\n{1}\n#modelling\n{2}\n#pipeline\n{3}\n'''.format(
            {key: value for key, value in content.items() if key in self.docking_form},
            {key: value for key, value in content.items() if key in self.scoring_function},
            {key: value for key, value in content.items() if key in self.modelling_form},
            content
            )
        # 加密
        content = self.make_secure(form, encryption=True)
        # 写文件
        with open(self.config_file, 'wb') as f:
            f.write(content)

    def submit_job(self, request, id_name, path_log):
        # 获取表单
        job_form = request.POST.dict()
        # 实例化model对象
        this_model = model_control(id_name)
        # 定义sh文件路径
        sh_file = '{}/pipeline.sh'.format(path_log)
        # 定义文件开头
        con = '''#!/bin/bash\nexport PYTHONPATH=$PYTHONPATH:/home/xujun/Project_1/mlplic\nmodule load anaconda3/5.1.0\n'''
        # 把表单记录到文件中传递参数
        del job_form['csrfmiddlewaretoken']  # 删除csrf键
        job_form['train_size'] = 0.8  # 设定超参数寻优时的训练集比例
        self.write_config(job_form)  # 写入文件
        # 判断是否对接
        if job_form['whether_dock'] != 'no':
            docked = False  # 赋予变量用于分割分子
            # 创建对接命令变量
            con += '''python {0} {1} {2} {3}\npython {4} {1}'''.format(self.split_py, id_name, self.job_type.strip('#'), docked,
                                                                       self.docking_py)
        else:
            docked = True
            # 创建对接命令变量
            con += '''python {0} {1} {2} {3}'''.format(self.split_py, id_name, self.job_type,
                                                       docked)
        # 循环判断被选择的描述符
        for sf in job_form:
            if sf in self.scoring_function:  # 若被选中
                # 创建执行命令变量
                new_con = '''\npython {}/{}.py {}'''.format(self.sfpy_path, self.scoring_function[sf], id_name)
                # 像sh中加入执行命令
                con += new_con
                # 创建对应文件夹
                os.mkdir('{}/{}/{}'.format(self.file_path, id_name, self.scoring_function[sf]))
        # 加入合并csv脚本
        con += '\npython {0} {1}\npython {0} {1} {2}'.format(self.merge_py, id_name, self.test_descriptors_csv)
        # 判断是否建模
        if job_form['algorithm'] != 'No Modelling':
            # 加入建模预测脚本
            con += '''\npython {}/{}.py {}'''.format(self.model_path, this_model.algorithm[job_form['algorithm']],
                                                     id_name)
            # 加入虚筛脚本
            con += '''\npython {} {}'''.format(self.screening_py, id_name)
        # 加入新的共同任务:压缩结果文件
        con = self.add_job(con, id_name, 'pipeline', request)
        # 创建sh文件
        with open(sh_file, 'w') as f:
            f.write(con)
        # 提交任务,并传入表单
        self.qsub(path_log, id_name, sh_file)
