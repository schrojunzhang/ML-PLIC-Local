from django.shortcuts import render
# Create your views here.
import os
from django.http import HttpResponse, HttpResponseRedirect
from base_scripts.base_class import global_base


def main(request):
    # 若提交表单，则进行计算
    if request.method == 'POST':
        # 生成任务id
        id_name = global_base.job_generation().get_id()
        # 实例化对象
        this_job = global_base.docking_control(id_name)
        # 创建任务路径
        path_job = this_job.path_job
        # 创建未对接文件存放路径
        path_files = this_job.path_for_dock
        # 创建日志路径
        path_log = this_job.path_for_log
        # 创建配置文件路径
        path_config = this_job.path_for_config
        # 创建报告路径
        path_report = this_job.path_for_report
        # 创建复合物路径
        path_complex = this_job.path_for_complex
        # 创建结果路径
        path_results = this_job.path_for_lig
        # 创建文件夹
        for dir in [path_job, path_files, path_log, path_results, path_config, path_report, path_complex]:
            os.mkdir(dir)
        # 循环上传csv小分子文件, 在前端要求上传文件, descriptor_file 必须上传，label_file可以没有
        for input_name in ['protein_file', 'ligand_file', 'crystal_file']:
            this_job.handle_uploaded_file(request, input_name, path_files)
        # 提交任务
        this_job.submit_job(request, id_name, path_log)
        # 重定向获取密钥
        return HttpResponseRedirect('../temp/{}'.format(id_name))
    # 未提交表单，渲染正常页面
    else:
        # 获取访问量
        this_job = global_base.job_generation()
        total_viewed = this_job.get_total_viewed()
        context = {'total_viewed': total_viewed}
        return HttpResponse(render(request, 'docking/index.html', context))
