#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 14:34:54 2021

@author: tristanpwdennis
"""

import allel; print('scikit-allel', allel.__version__)
import sys
import numcodecs
import argparse

#args and pod
parser = argparse.ArgumentParser(prog='vcf_to_zarr',description='convert a vcf file to ZARR format')
parser = argparse.ArgumentParser()

parser.add_argument('-v', '--vcf', type=str, required=True, help='the number of samples you want to generate')
parser.add_argument('-z', '--zarr_prefix', type=str, required=True, help='prefix of your outfile: prefix.zarr')

args = parser.parse_args()

allel.vcf_to_zarr(args.vcf, args.zarr_prefix,  fields = '*', alt_number=8, log=sys.stdout, compressor=numcodecs.Blosc(cname='zstd', clevel=1, shuffle=False),  exclude_fields='variants/numalt', overwrite=True) 


