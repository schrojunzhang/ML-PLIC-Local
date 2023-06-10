#!usr/bin/env python
# -*- coding:utf-8 -*-
from django.urls import path
from . import views
app_name = 'cal_model'
urlpatterns = [
    path('modelling/', views.main, name='modelling')
]