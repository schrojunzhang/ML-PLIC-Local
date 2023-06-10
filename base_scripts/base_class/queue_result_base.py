#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-18

import os
import glob
import shutil
import time
import pandas as pd
import numpy as np
from multiprocessing import Pool
from base_scripts.base_class import global_base
from django.core.paginator import Paginator
from django.http import FileResponse


# 队列系统
class queue(global_base.job_generation):

    def __init__(self):
        super().__init__()
        self.rsa_file = '{}_rsa_file'  # 密钥文件
        # 任务类型分类
        self.job2time = {'#docking': self.docking_progress,
                         '#descriptors': self.descriptors_progress,
                         '#modelling': self.modelling_progress}
        # 下载文件与文件后缀的映射
        self.download_file_extension = {
            'protein': '.pdb',
            'ligand': '.mol2',
            'decoys': '.mol2',
            'crystal_ligand': '.mol2',
            'test': '.mol2',
            'descriptor': '.zip',
            'screening': '.zip',
            'descriptor_file': '.csv',
            'test_descriptor_file': '.csv',
            'config': '',
        }
        # 密码到结果文件路径的转义
        self.str2code = {
            '+': '%2B',
            ' ': '%20',
            '/': '%2F',
            '?': '%3F',
            '%': '%25',
            '#': '%23',
            '&': '%26',
            '=': '%3D'
        }
        # example 文件下载
        self.path_for_download = '/home/xujun/Project_1/mlplic/files_for_download'
        # 待下载文件
        self.file_for_download = self.path_for_download + '/{}{}'

    def get_jobs(self):  # 获取所有任务
        return glob.glob('{}/*/log/*.sh'.format(self.file_path))

    # 进度条
    def cal_percent(self, total, now_value):
        try:  # 防止出现被除数为0的情况 比如在pipline中，不选择对接时，但是config中有#docking标志，会导致除0
            percent = now_value / total
        except:  # 出错则直接返回1
            percent = 1
        return percent

    def docking_progress(self, job_id):
        # 实例化对象  根据分子数量进行进度判断
        this_job = global_base.docking_control(job_id)
        # 判断是否已经对接完毕
        if os.path.exists(this_job.path_for_dock):  # 未对接完毕
            # 判断准备进度
            total = len(glob.glob('{}/*.mol2'.format(this_job.path_for_dock)))  # 未准备小分子总数  mol2的数量
            now_value = len((glob.glob('{}/*.pdbqt'.format(this_job.path_for_dock)))) - 1  # 准备后的小分子数量 pdbqt的数量
            prep_progress = self.cal_percent(total, now_value)
            # 判断对接进度, 文件大小大于1KB
            now_value = len(
                [1 for file in glob.glob('{}/*.mol2'.format(this_job.path_for_lig)) if os.path.getsize(file) > 1000])
            # 对接进度
            docking_progress = max(0.2 * prep_progress + 0.8 * self.cal_percent(total, now_value), 0.01)
        else:  # 已对接完毕,直接返回1
            docking_progress = 1
        # 返回数据
        return docking_progress

    def descriptors_progress(self, job_id):
        # 根据文件夹数量进行进度判断
        # 实例化对象
        this_job = global_base.job_control(job_id)
        # 待计算描述符数量
        job_form = this_job.read_config(this_job.job_type)  # 获取描述符表单
        total = len([1 for sf in this_job.scoring_function if sf in job_form])  # 获取应该计算的描述符数量
        wait_value = len([1 for sf in this_job.scoring_function if
                          os.path.exists(
                              '{}/{}'.format(this_job.path_job, this_job.scoring_function[sf]))])  # 获取还未计算的描述符
        now_value = total - wait_value  # 获取已经计算的描述符
        # 返回
        return self.cal_percent(total, now_value)

    def modelling_progress(self, job_id):
        # 实例化模型对象
        this_model = global_base.model_control(job_id)
        # 获取表单
        job_form = this_model.read_config(job_type=this_model.job_type)
        # 获得寻优次数
        max_eval = job_form.get('max_eval')
        # 判断pipeline是否选择建模
        if max_eval:  # 选择了建模, 转成整数
            max_eval = int(max_eval)
        else:  # 未选择建模，直接返回进度
            return 0.95
        # 若寻优次数大于0
        if max_eval > 0:
            # 判断是否生成了CSV文件
            if not os.path.exists(this_model.bias_var_csv):  # 未生成CSV，说明还在前面的数据读取和预处理中
                # 进度赋予0.01
                progress = 0.01
            else:  # 生成了CSV，开始寻优
                # 读取方差偏差分析数据csv
                df = pd.read_csv(this_model.bias_var_csv, encoding='utf-8', header=None)
                # 获取已经寻优的次数
                now_eval = df.shape[0]
                # 计算进度
                progress = 0.3 + 0.7 * now_eval / max_eval - 0.05
            # 返回进度
            return progress
        else:  # 若寻优次数为0
            return 0.95

    def detail_status(self, job_id):  # 获取进度条
        # 判断任务类型
        this_job = global_base.job_control(job_id)  # 实例化任务对象
        # 初始化总进度
        total_progress = 0
        # 初始化存在的任务数
        job_num = 0
        # 循环判断有哪些任务类型
        for job_ in self.job2time:
            # 读取表单查看是否存在该任务
            if this_job.read_config(job_):  # 若存在任务
                total_progress += self.job2time[job_](job_id)  # 加上进度条
                job_num += 1  # 任务数+1
        # 进度条
        total_progress = total_progress / job_num * 100 - 1  # -1 任务进度条按照任务数均摊，并放大到100，最后-10， 保证进度条最大为99
        # 控制进度条最小值为0
        if total_progress <= 0:
            total_progress = 1
        # 返回进度条
        return total_progress

    # 获取单个任务信息
    def get_job_info(self, job_sh):  # 获取单个任务的信息
        # 获取任务名
        job_id = job_sh.split('/')[-3]
        # 实例化单个任务对象
        this_job = global_base.job_control(job_id)
        # 单个任务的日志路径
        path_log = this_job.path_for_log
        # 获取任务类型
        # sh_file = glob.glob('{}/*.sh'.format(path_log))  # 获取sh文件
        # # 判断文件是否存在
        # if len(sh_file) == 0:
        #     # 文件不存在，任务提交失败
        #     return None
        # else:
        #     job_type = sh_file[0].split('/')[-1].split('.')[0]  # 获取任务类型
        job_type = job_sh.split('/')[-1].split('.')[0]  # 获取任务类型
        # 获取任务时间
        time_number = os.path.getctime(job_sh)  # 时间戳
        job_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time_number))  # 有格式的时间
        # 获取任务状态
        error_file = glob.glob('{}/*.e*'.format(path_log))  # 错误文件
        submit_time_log = '{}/{}_time.log'.format(path_log, job_id)  # 时间log
        result_file = this_job.result_file  # 结果文件
        # 判断状态
        if os.path.exists(result_file):  # 若结果文件存在，则任务完成
            job_status = 100
        elif len(error_file) == 1:  # 若结果文件不存在，但是报错文件存在，则任务出错
            job_status = -1
        elif os.path.exists(submit_time_log) == 1:  # 若时间文件存在，则任务开始计算
            job_status = self.detail_status(job_id)
        else:  # 若上述文件都不存在，则任务还在排队
            job_status = 0
        # 给出链接
        job_result_url = self.make_secure(content=job_id, encryption=True)  # 给出加密后的链接
        # 整合结果
        job_info = [job_id, job_type, job_time, job_status, job_result_url, time_number]
        # 增加结果至列表中
        return job_info

    # 分页
    def split_pages(self, request, job_infos):
        # 实例化分页器
        asfp_paginator = Paginator(job_infos, 15)
        # 从表单获取页数
        page_num = request.GET.get('page')
        # 生成对应页数的数据, 若输入的页数格式不对，则直接返回第一页数据
        try:
            page_obj = asfp_paginator.get_page(page_num)
        except:
            page_obj = asfp_paginator.get_page('1')
        # 返回数据
        return page_obj

    # 创建密钥
    def write_key(self, job_id):
        # 加密
        content = self.make_secure(job_id, encryption=True)
        # 写文件
        with open(self.key_file.format(job_id), 'wb') as f:
            f.write(content)

    def download_file(self, dst_file, file_name):
        file = open(dst_file, 'rb')  # 打开文件
        response = FileResponse(file)  # 实例化对象
        response['Content-Type'] = 'application/octet-stream'  # 定义文件传输类型
        response['Content-Disposition'] = 'attachment;filename={}'.format(file_name)  # 定义文件组成
        return response

    # 获取结果时验证密码
    def check_passwd(self, request, job_id):
        # 获取上传的文件
        rsa_file = request.FILES[self.rsa_file.format(job_id)]
        # 判断文件格式是否正确，只有.rsd可以被阅读
        file_format = rsa_file.name.split('.')[-1]
        # 若文件格式正确
        if file_format == 'rsa':
            # 读取文件数据
            content = rsa_file.read()
            # 解密
            rsa_id = self.make_secure(content, encryption=False)
            # 验证密码
            if job_id == rsa_id:
                return True
            else:
                return False
        else:  # 文件格式错误
            return False

    # 每隔1h检查是否任务是否完成， 先判断是否上传了邮箱，先判断result.zip是否存在，若pbs文件也存在，则发送邮件，发送后删除pbs文件，
    # 若不存在，则不发送
    def finished_jobs(self, job_id):
        # 实例化单个任务对象
        this_job = global_base.job_control(job_id)
        # 判断email是否存在
        job_form = this_job.read_config('#')  # 读取第一个表单
        email = job_form.get('my_email', None)  # 从表单获取email
        if email:  # 若存在email
            # 获取任务状态
            result_file = this_job.result_file  # 结果文件
            pbs_file = this_job.pbs_file  # pbs文件
            # 任务若是完成
            if os.path.exists(result_file) and os.path.exists(pbs_file):  # 若两个文件同时存在
                this_job.send_email(email, job_id, this_job.job_finished_content,
                                    ath_file=this_job.key_file.format(job_id))  # 发送邮件
                os.remove(pbs_file)  # 删除pbs文件，防止下次继续发送邮件


# 处理结果类
class result(global_base.job_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        # 结果路径
        self.path_for_result = '{}/result'.format(self.path_job)
        # 根据任务类型不同，生成不同结果文件
        self.job2result = {
            'docking': [self.path_for_lig, self.path_for_report, self.path_for_complex],
            'descriptors': [self.path_for_csv],
            'modelling': [self.path_for_report, self.path_for_config],
            'screening': [self.path_for_report],
            'pipeline': [self.path_for_lig, self.path_for_csv, self.path_for_complex, self.path_for_report]
        }
        # 根据不同任务类型，删除不同文件夹
        self.job2del = {
            'docking': [self.path_for_lig, self.path_for_result, self.path_for_dock],
            'descriptors': [self.path_for_csv, self.path_for_report, self.path_for_lig, self.path_for_result],
            'modelling': [self.path_for_csv, self.path_for_result],
            'screening': [self.path_for_csv, self.path_for_result],
            'pipeline': [self.path_for_lig, self.path_for_csv, self.path_for_result]
        }
        # score.csv文件中的相互作用前缀
        self.interaction_types = {'hb': 'hydrogen bond',
                                  'qq': 'salt bridge',
                                  'lipo': 'hydrophobic contacts',
                                  'clb': 'halogenbond',
                                  'metal': 'metal'}
        # deal_result.py 中根据任务类型不同，读取不同的小分子名字CSV
        self.type2ligandcsv = {
            'pipeline': self.active_ligand_csv,
            'docking': self.docking_score_csv
        }
        # 用于生成最大子结构的mol 原始小分子的路径
        self.ligs_for_mcs_src = self.path_for_lig + '/{}.mol2'
        # 对接任务完成后合并成一个配体
        self.ligs_docked = '{}/docked_ligands.mol2'.format(self.path_for_report)
        # 生成最大子结构的mol小分子的路径
        self.ligs_for_mcs_dst = self.path_for_dock + '/{}.mol'
        # 聚类脚本路径
        self.cluster_py = '{}/ligand_cluster.py'.format(self.help_path)
        # 聚类结构图
        self.cluster_png = '{}/cluster.png'.format(self.path_for_report)
        # 用于计算相互作用的脚本
        self.interaction_py = '{}/cal_interactions.py'.format(self.help_path)

    # 获取任务类型
    def get_job_type(self):
        # 单个任务的日志路径
        path_log = self.path_for_log
        # 获取任务类型
        sh_file = glob.glob('{}/*.sh'.format(path_log))  # 获取sh文件
        job_type = sh_file[0].split('/')[-1].split('.')[0]  # 获取任务类型
        return job_type

    # 获取得分前100的小分子
    def top100_ligs(self, csv):
        # 读取文件
        df = pd.read_csv(csv, encoding='utf-8')
        # 排序
        if df.score.sum() > 0:
            df.sort_values(by='score', ascending=False, inplace=True)
        else:
            df.sort_values(by='score', inplace=True, ascending=True)
        # 取前100的小分子
        lig_names = df.iloc[:100, 0].tolist() 
        # 返回前100的分子
        return lig_names

    def cal_cluster(self, active_names):
        # 通过描述符文件是否存在判断是对接任务还是pipeline任务
        if os.path.exists(self.descriptor_csv):  # pipeline
            df = pd.read_csv(self.descriptor_csv, encoding='utf-8')
            # 获取活性分子
            lig_names = df[df.iloc[:, -1] == 1].iloc[:, 0].values
            active_names.extend(lig_names)
        # 转换分子格式,用于计算最大子结构
        src_dst_ligs = zip([self.ligs_for_mcs_src.format(i) for i in active_names],
                           [self.ligs_for_mcs_dst.format(i) for i in active_names])
        pool = Pool(28)
        pool.starmap(self.openbabel_transform, src_dst_ligs)
        pool.close()
        pool.join()
        # 计算最大子结构
        cmd = '{} {} {} {} {} '.format(self.ifp_module, self.cluster_py, self.path_for_complex, self.path_for_dock,
                                       self.cluster_png)
        os.system(cmd)

    def cal_interactions(self, job_id):
        # cmd
        cmd = '{} {} {}'.format(self.ifp_module, self.interaction_py, job_id)
        # 执行
        os.system(cmd)

    def read_data_csv(self, prefix):
        # 按行读取
        with open(self.data_collect_csv, 'r') as f:
            con = f.readlines()
        # 循环判断
        for i in range(len(con)):
            if con[i].startswith(prefix):
                data = con[i + 1]  # 获取下一行数据列表
                data = eval(data)  # 转成列表
                # 判断是否有数据
                if data:
                    data = [prefix, data]  # 首位插入作用名字
                    return data
        # 若没有找到数据
        return None

    def merge_lig(self):
        # 初始化内容
        content = ''''''
        # 获取对接后的小分子
        mol2s = glob.glob('{}/*.mol2'.format(self.path_for_lig))
        # 合并分子
        for mol2 in mol2s:
            with open(mol2, 'r') as f:
                content += '{}\n'.format(f.read())
        # 创建分子
        with open(self.ligs_docked, 'w') as f:
            f.write(content)

    # 根据任务类型不同生成不同结果文件
    def zip_result(self):
        # 创建result 文件夹
        os.mkdir(self.path_for_result)
        # 获取任务类型
        job_type = self.get_job_type()
        # 移动文件夹到result文件夹内
        for dir_ in self.job2result[job_type]:
            # 命令行
            cmd = 'cp -r {} {}'.format(dir_, self.path_for_result)
            # 执行
            os.system(cmd)
        # 压缩文件夹
        cmd = 'cd {}&&zip -r result.zip result'.format(self.path_job)
        # 执行压缩
        os.system(cmd)
        # 删除无用文件夹
        for dir_ in self.job2del[job_type]:
            # shutil.rmtree(dir_)
            pass

    # 处理对接
    def deal_docking(self, csv=None):
        # 若没有给出CSV
        if not csv:
            csv = self.docking_score_csv
            # 收集前100的分数  docking
            df = pd.read_csv(csv).sort_values(by='score')  # 读取csv
        else:
            # 收集前100的分数  pipline
            df = pd.read_csv(csv).sort_values(by='score', ascending=False)  # 读取csv
        # 取100， 不足一百则全取完
        lig_infos = zip([i for i in df.iloc[:100, 0].tolist()], [j for j in df.iloc[:100, 1].tolist()])
        # 获取每个分子的相互作用信息
        df_each = pd.read_csv(self.interaction_csv, encoding='utf-8').dropna()
        each_interaction_static = {}
        # 循环获取每个分子相互作用信息
        for i in range(df_each.shape[0]):
            interactions = df_each.iloc[i, 1:].values.tolist()
            # 合并成字典
            each_interaction_static[df_each.iloc[i, 0]] = [eval(j) for j in interactions]
        # 获取所有活性分子相互作用统计信息
        interaction_res_frequency = []  # 初始化空列表
        for interaction_type in self.interaction_types:
            # 读取数据
            df = pd.read_csv(self.data_collect_csv, encoding='utf-8', header=None)
            # 读取数据后为numpy.array类, 直接转成列表
            rec_frequnece = df[df.iloc[:, 0] == interaction_type].iloc[:, 1].values.tolist()
            # 若数据不为空, 过滤不存在的相互作用类型
            if rec_frequnece != [np.nan]:
                # 把二重列表字符串转换成列表
                rec_frequnece = [eval(i) for i in rec_frequnece]
                # 给列表加上作用名称换名字
                rec_frequnece.insert(0, self.interaction_types[interaction_type])  # 换名字
                # 增加列表
                interaction_res_frequency.append(rec_frequnece)
        # 返回信息 smarts,
        return lig_infos, each_interaction_static, interaction_res_frequency

    # 处理描述符任务 返回描述符文件
    def result_download(self):
        file = open(self.result_file, 'rb')  # 打开文件
        response = FileResponse(file)  # 实例化对象
        response['Content-Type'] = 'application/octet-stream'  # 定义文件传输类型
        response['Content-Disposition'] = 'attachment;filename={}'.format('result.zip')  # 定义文件组成
        return response

    # 处理建模
    def deal_modelling(self):
        # 收集评价指标
        model_trics = pd.read_csv(self.model_metric_csv, encoding='utf-8').iloc[0, :].tolist()
        # 收集特征重要性
        df = pd.read_csv(self.feature_importance_csv, encoding='utf-8')
        features = df.columns.tolist()
        importance = df.iloc[0, :].tolist()
        # 最终参数
        df = pd.read_csv(self.final_parameter_csv, encoding='utf-8')
        parameter_names = df.columns.tolist()
        parameters = df.iloc[0, :].tolist()
        # 返回参数
        return model_trics, features, importance, parameter_names, parameters

    # 处理虚筛
    def deal_screening(self):
        self.result_download()

    # 处理pipeline
    def deal_pipline(self):
        # 获取小分子数据 smarts,
        lig_infos, each_interaction_static, interaction_res_frequency = self.deal_docking(csv=self.active_ligand_csv)
        # 获取模型数据
        model_trics, features, importance, parameter_names, parameters = self.deal_modelling()
        # 返回参数
        return model_trics, features, importance, parameter_names, parameters, lig_infos, each_interaction_static, interaction_res_frequency
