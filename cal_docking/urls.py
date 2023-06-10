#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-15

from django.urls import path
from . import views
app_name = 'cal_docking'
urlpatterns = [
    path('docking/', views.main, name='docking')
]