#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-18
from django.urls import path, re_path
from . import views
app_name = 'cal_queue'
urlpatterns = [
    path('queue', views.queue_page, name='queue'),  # 队列
    path('help', views.help, name='help'),  # 帮助
    path('paper', views.paper, name='paper'),  # 文章
    path('contact', views.contact, name='contact'),  # 联系
    path('viewed', views.viewed, name='viewed'),  # IP访问次数
    re_path(r'^results/(?P<job_id>job_[0-9]{8})', views.check_key, name='result'),  # 验证表单
    re_path(r'^results/(?P<job_id>b.*cadd$)', views.show_result, name='result_show'),  # 结果界面
    re_path(r'^results/(?P<job_id>b.*cadd)/download', views.download_result, name='download'),  # 结果界面下载结果文件
    re_path(r'^temp/job.*', views.temp, name='temp'),  # 各个任务完成后转跳的空白界面
    re_path(r'^temp/download/(?P<job_id>.*)', views.generate_key, name='generate_key'),  # 空白界面下载密钥
    re_path(r'^example/(?P<file_name>.*)', views.example_download, name='example_download'),  # 下载example 文件
]