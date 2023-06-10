from django.shortcuts import render

# Create your views here.
import os
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect
from base_scripts.base_class import global_base, scoring_function_base, fingerprint_base


def main(request):
    # 若提交表单，则进行计算
    if request.method == 'POST':
        # 生成任务id
        id_name = global_base.job_generation().get_id()
        # 实例化对象
        this_job = global_base.job_control(id_name)
        # 创建任务路径
        path_job = this_job.path_job
        # 创建分子文件存放路径
        path_files = this_job.path_for_lig
        # 创建日志路径
        path_log = this_job.path_for_log
        # 配置文件路径
        path_config = this_job.path_for_config
        # 创建报告路径
        path_report = this_job.path_for_report
        # 创建结果路径
        path_results = this_job.path_for_csv
        # 创建文件夹
        for dir in [path_job, path_files, path_log, path_results, path_config, path_report]:
            os.mkdir(dir)
        # 循环上传蛋白小分子文件, 在前端要求上传文件，这里就不再检查是否包含文件[晶体结构还需检查]
        for input_name in ['protein_file', 'ligand_file', 'crystal_file', 'decoys_file', 'test_file']:
            this_job.handle_uploaded_file(request, input_name, path_files)
        # 提交任务
        this_job.submit_job(request, id_name, path_log)
        # 重定向获取密钥
        return HttpResponseRedirect('../temp/{}'.format(id_name))

    # 未提交表单，渲染正常页面
    else:
        # 实例化sf对象
        this_sf = scoring_function_base.scoring_function('')
        # 获取访问量
        total_viewed = this_sf.get_total_viewed()
        # 实例化指纹对象
        this_fp = fingerprint_base.fingerpringt('')
        this_ifp = fingerprint_base.interaction_fingerprint('')
        # 获取字典并化简
        sf_dic = this_sf.sf2energy_term
        # 化简sf
        sf_dic['nn_vina'] = ['nn_348']
        sf_dic['rfscore_sybyl'] = ['sybyl_5046']
        sf_dic['rfscore_element'] = ['element_486']
        # 化简fp
        fp_dic = {i: ['{}_{}'.format(i, len(j)-1)] for i, j in this_fp.fp2fingerprint.items()}
        ifp_dic = {i: ['{}_{}'.format(i, len(j)-1)] for i, j in this_ifp.ifp2fingerprint.items()}
        # 合并对象
        sf_descriptors = dict(**sf_dic, **fp_dic, **ifp_dic)
        # 获取描述符数据
        context = {'sfs': sf_descriptors.items(), 'total_viewed': total_viewed}
        return render(request, 'descriptors/index.html', context)
