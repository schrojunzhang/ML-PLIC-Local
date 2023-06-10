#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-16

import sys
from base_scripts.base_class.global_base import job_control

def main():
    # 获得任务ID
    job_id = sys.argv[1]
    # 实例化对象
    this_job = job_control(job_id)
    # 获得目标CSV
    try:
        dst_file = sys.argv[2]
    except:
        dst_file = this_job.descriptor_csv
    # 合并csv
    this_job.merge_csv(dst_file)

if __name__ == '__main__':
    main()
