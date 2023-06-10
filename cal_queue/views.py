from django.shortcuts import render

# Create your views here.
import os
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect, Http404
from base_scripts.base_class import queue_result_base, global_base
from django.contrib import messages
from multiprocessing import Pool
import re


def temp(request):
    # 返回页面
    return render(request, 'temp/index.html')


def viewed(request):
    # 重定向
    return HttpResponseRedirect('https://www.revolvermaps.com/livestats/5z9d6s1pvmv/')


def paper(request):
    # 获取访问量
    this_job = queue_result_base.queue()
    total_viewed = this_job.get_total_viewed()
    context = {'total_viewed': total_viewed}
    # 返回页面
    return render(request, 'paper/index.html', context)


def help(request):
    # 获取访问量
    this_job = queue_result_base.queue()
    total_viewed = this_job.get_total_viewed()
    context = {'total_viewed': total_viewed}
    # 返回页面
    return render(request, 'help/index.html', context)


def contact(request):
    # 获取访问量
    this_job = queue_result_base.queue()
    total_viewed = this_job.get_total_viewed()
    context = {'total_viewed': total_viewed}
    # 返回页面
    return render(request, 'contact/index.html', context)


# 在temp页面生成密钥并返回下载，若存在邮箱，则把文件发送给邮箱
def generate_key(request, job_id):
    # 实例化对象
    this_result = queue_result_base.queue()
    # 获取密钥
    key_file = this_result.key_file.format(job_id)
    # 创建密钥
    this_result.write_key(job_id)
    # 实例化任务控制对象
    this_job = global_base.job_control(job_id)
    # 读取第一个表单
    job_form = this_job.read_config(job_type='#')
    # 若邮箱存在，则发送邮件
    if job_form.get('my_email', None):
        this_job.send_email(job_form['my_email'], job_id, this_job.email_submit_content, key_file)
    # 返回文件
    return this_result.download_file(dst_file=key_file, file_name='{}.rsa'.format(job_id))


# 每隔1h检测并判断是否发送邮件
def check_send_email():
    # 实例化对象
    this_queue = queue_result_base.queue()
    # 获取任务名
    job_shs = this_queue.get_jobs()
    job_ids = [i.split('/')[-3] for i in job_shs]
    # 多进程任务信息
    pool = Pool(28)
    pool.map(this_queue.finished_jobs, job_ids)
    pool.close()
    pool.join()


# 根据上传的密钥验证
def check_key(request, job_id):
    # 若提交表单，则进行计算
    if request.method == 'POST':
        # 实例化对象
        this_result = queue_result_base.queue()
        # 对比密码
        if this_result.check_passwd(request, job_id=job_id):  # 若密码正确
            # 对任务名加密形成结果网页的URL
            job_result_url = str(this_result.make_secure(content=job_id, encryption=True))
            # 循环替换转义字符
            for string, code in this_result.str2code.items():
                job_result_url = job_result_url.replace(string, code)
            # 判断末尾加上cadd，作为参数截至的标志
            job_result_url = job_result_url + 'cadd'
            # 重定向到结果页面
            return HttpResponseRedirect(job_result_url)
        else:  # 若密码不正确
            # 弹窗
            messages.success(request, 'Validation Failure!')
            # 重定向回队列页面
            return HttpResponseRedirect('../queue')
    # 未提交表单，渲染空白页面
    else:
        return HttpResponseRedirect('../queue')


#  获取对应链接后渲染对应结果界面
def show_result(request, job_id):
    # 实例化对象
    this_result = queue_result_base.queue()
    # 换回问号
    # 循环把转义字符换回正常符号
    for string, code in this_result.str2code.items():
        job_id = job_id.replace(code, string)
    job_id = job_id.replace('cadd', '')  # 去掉后缀
    # 解密
    job_id = eval(job_id)  # 从字符串转成二进制
    job_id = this_result.make_secure(job_id, encryption=False)  # 解密
    # 获取任务类型
    this_result = queue_result_base.result(job_id)
    job_type = this_result.get_job_type()
    # 获取访问量
    total_viewed = this_result.get_total_viewed()
    # 先判断是否为pipeline
    if job_type == 'pipeline':
        # 实例化pipeline对象
        this_pipeline = global_base.pipline_control(job_id)
        # 判断是否对接和建模
        docking = this_pipeline.whether_dock()
        modelling = this_pipeline.whether_modelling()
        # 建模
        if modelling:
            model_metrics, feature, importance, parameter_names, parameters, lig_infos, each_interaction_static, interaction_res_frequency = this_result.deal_pipline()  # 获取信息
            # 获取数据
            context = {'job_type': job_type, 'job_id': job_id, 'lig_infos': lig_infos,
                       'each_interaction_static': each_interaction_static,
                       'interaction_rec_frequency': interaction_res_frequency, 'model_metrics': model_metrics,
                       'feature': feature, 'importance': importance, 'parameter_names': parameter_names,
                       'parameters': parameters, 'total_viewed': total_viewed}
            return render(request, 'result/index.html', context)  # 对应的路径为 /asfp/results/... 会报错，但是不影响渲染
        elif docking:  # 若不建模但是对接
            job_type = 'docking'  # 任务类型改为对接
        else:  # 不建模且不对接
            job_type = 'descriptors'  # 任务类型改为描述符提取
    # 若pipline未建模或者不是pipeline，再判断其他任务
    if job_type == 'docking':  # 对接 或者 pipeline进行对接未建模
        lig_infos, each_interaction_static, interaction_res_frequency = this_result.deal_docking()  # 获取信息
        # 获取数据
        context = {'job_type': job_type, 'job_id': job_id, 'lig_infos': lig_infos,
                   'each_interaction_static': each_interaction_static,
                   'interaction_rec_frequency': interaction_res_frequency, 'total_viewed': total_viewed}
        return render(request, 'result/index.html', context)  # 对应的路径为 /asfp/results/... 会报错，但是不影响渲染
    elif job_type == 'descriptors':  # 计算描述符 返回结果文件  descriptor 或者 pipeline未对接未建模
        return this_result.result_download()
    elif job_type == 'modelling':  # 建模
        model_metrics, feature, importance, parameter_names, parameters = this_result.deal_modelling()  # 获取信息
        # 获取数据
        context = {'job_type': job_type, 'job_id': job_id, 'model_metrics': model_metrics, 'feature': feature,
                   'importance': importance, 'total_viewed': total_viewed,
                   'parameter_names': parameter_names, 'parameters': parameters}
        return render(request, 'result/index.html', context)  # 对应的路径为 /asfp/results/... 会报错，但是不影响渲染
    elif job_type == 'screening':  # 虚筛 返回结果文件
        return this_result.result_download()


# 结果界面获取下载文件
#  获取对应链接后渲染对应结果界面
def download_result(request, job_id):
    # 实例化对象
    this_result = queue_result_base.queue()
    # 循环把转义字符换回正常符号
    for string, code in this_result.str2code.items():
        job_id = job_id.replace(code, string)
    job_id = job_id.replace('cadd', '')  # 去掉后缀
    # 解密
    job_id = eval(job_id)  # 从字符串转成二进制
    job_id = this_result.make_secure(job_id, encryption=False)  # 解密
    # 实例化任务对象
    this_result = queue_result_base.result(job_id)
    # 下载文件
    return this_result.result_download()


# 排队系统
def queue(request):  # 获取渲染页面所需数据
    # 实例化对象
    this_queue = queue_result_base.queue()
    # 获取访问量
    total_viewed = this_queue.get_total_viewed()
    # 根据sh文件获取任务名
    job_shs = this_queue.get_jobs()
    # 多进程任务信息
    pool = Pool(28)
    job_infos = pool.map(this_queue.get_job_info, job_shs)
    pool.close()
    pool.join()
    # 根据时间信息排序
    job_infos.sort(key=lambda x: x[-1], reverse=True)
    # 获取分页信息
    page_obj = this_queue.split_pages(request, job_infos)
    # 返回给页面的信息
    context = {'jobs': page_obj.object_list, 'page_obj': page_obj, 'start': page_obj.start_index() - 1,
               'total_viewed': total_viewed}
    # 返回
    return context


def queue_page(request):  # 渲染排队页面
    # 获取数据
    context = queue(request)
    # 渲染页面
    return render(request, 'queue/index.html', context)


# example 文件下载
def example_download(request, file_name):
    # 实例化对象
    this_queue = queue_result_base.queue()
    # 检查文件名是否正确
    if file_name not in this_queue.download_file_extension:
        raise Http404('Bad URL')
    # 定义需要下载的文件
    dst_file = this_queue.file_for_download.format(file_name, this_queue.download_file_extension[file_name])
    # 文件名字
    file_name = 'example_{}'.format(dst_file.split('/')[-1])
    # 下载文件
    return this_queue.download_file(dst_file=dst_file, file_name=file_name)
