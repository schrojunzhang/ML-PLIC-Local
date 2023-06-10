from django.shortcuts import render

# Create your views here.
import os
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect
from base_scripts.base_class import global_base
from base_scripts.base_class import scoring_function_base, fingerprint_base


def main(request):
    # 若提交表单，则进行计算
    if request.method == 'POST':
        # 生成任务id
        id_name = global_base.job_generation().get_id()
        # 实例化对象
        this_job = global_base.pipline_control(id_name)
        # 创建任务路径
        path_job = this_job.path_job
        # 创建对接路径
        path_undocked = this_job.path_for_dock
        # 创建文件路径
        path_files = this_job.path_for_lig
        # 创建csv文件存放路径
        path_csvs = this_job.path_for_csv
        # 创建复合物路径
        path_complex = this_job.path_for_complex
        # 创建日志路径
        path_log = this_job.path_for_log
        # 创建报告路径
        path_report = this_job.path_for_report
        # 创建配置文件路径
        path_config = this_job.path_for_config
        # 创建文件夹
        for dir in [path_job, path_undocked, path_files, path_csvs, path_complex, path_log, path_config, path_report]:
            os.mkdir(dir)
        # 判断是否需要进行对接
        if request.POST.get('whether_dock') != 'no':
            file_path = path_undocked  # 判断文件上传路径
        else:
            file_path = path_files
        # 循环上传csv小分子文件, 在前端要求上传文件, 'protein_file', 'ligand_file'必须上传，其他可以没有
        for input_name in ['protein_file', 'ligand_file', 'crystal_file', 'decoys_file', 'test_file']:
            this_job.handle_uploaded_file(request, input_name, file_path)
        # 提交任务
        this_job.submit_job(request, id_name, path_log)
        # 重定向获取密钥
        return HttpResponseRedirect('temp/{}'.format(id_name))
    # 未提交表单，渲染正常页面
    else:
        # 实例化sf对象
        this_sf = scoring_function_base.scoring_function('')
        # 获取访问次数
        total_viewed = this_sf.get_total_viewed() + 1
        with open(this_sf.total_viewed, 'w') as f:
            f.write(str(total_viewed))
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
        fp_dic = {i: ['{}_{}'.format(i, len(j) - 1)] for i, j in this_fp.fp2fingerprint.items()}
        ifp_dic = {i: ['{}_{}'.format(i, len(j) - 1)] for i, j in this_ifp.ifp2fingerprint.items()}
        # 合并对象
        sf_descriptors = dict(**sf_dic, **fp_dic, **ifp_dic)
        # 获取描述符数据
        context = {'sfs': sf_descriptors.items(), 'total_viewed': total_viewed}
        return HttpResponse(render(request, 'pipline/index.html', context))
