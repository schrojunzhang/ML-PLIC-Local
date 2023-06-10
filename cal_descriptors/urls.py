#!usr/bin/env python
# -*- coding:utf-8 -*-
from django.urls import path
from . import views
app_name = 'cal_descriptors'
urlpatterns = [
    path('descriptors/', views.main, name='descriptors')
]