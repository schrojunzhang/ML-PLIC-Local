#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-11-20
# 用于分割分子
import sys
from base_scripts.base_class.global_base import job_control, docking_control, pipline_control

# job_id
job_id = sys.argv[1]
# job_type
job_type = sys.argv[2]
# docking
if job_type == 'docking':
    # instance
    this_job = docking_control(job_id)
    # split
    this_job.split_file(src_ligand=this_job.ligand, dst_path=this_job.path_for_dock, return_name=False)
    # descriptor
elif job_type == 'descriptors':
    # instance
    this_job = job_control(job_id)
    # split
    this_job.split_addLabel()
    # pipeline
else:
    # instance
    this_job = pipline_control(job_id)
    # get data
    docked = sys.argv[3]
    # split
    this_job.split_files(docked=docked)
