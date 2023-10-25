#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:36:26 2023

@author: rcruces
"""

import os
import glob
import nibabel as nib

subj_func='/data_/mica3/BIDS_PNI/derivatives/micapipe_v0.2.0/'

gii = glob.glob(subj_func+'sub-*/ses-*/func/desc-me_task-*/volumetric/*SNR*.shape.gii')

for gii_file in gii:
    # Split the string by '/'
    parts = gii_file.split('/')[11].replace('_space-func_desc-me_tSNR.shape.gii','')
    task = gii_file.split('/')[9].replace('desc-me_task-','')
    data = nib.load(gii_file).darrays[0].data
    print(f"{data.shape} : {parts} {task}")