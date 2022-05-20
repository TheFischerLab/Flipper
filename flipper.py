#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    FLIPPER: A program for identifying changes in side chain conformations 
#    from Ringer measurements
#    Authors: Tim Stachowski & Marcus Fischer
#    Email: tim.stachowski@stjude.org
#    Copyright 2022 St. Jude Children's Research Hospital, Memphis, TN
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import pandas as pd
import argparse
import time
from ringer_tools import *

CLI = argparse.ArgumentParser()

CLI.add_argument(
    '-i1',
    '--filename1',
    nargs="?",
    type=str,
    default=None,
    help='input file #1 - must be a standard Ringer output CSV with only one chain'
)

CLI.add_argument(
    '-i2',
    '--filename2',
    nargs="?",
    type=str,
    default=None,
    help='input file #2 - must be a standard Ringer output CSV with only one chain'
)

CLI.add_argument(
    '-t',
    '--sigmathreshold',
    nargs=1,
    type=float,
    default=0.3,
    help='sigma threshold for peak finding. default = 0.3'
)

CLI.add_argument(
    '-plot',
    '--plot',
    nargs=1,
    type=str,
    default=False,
    help='Save individual plots showing peak finding results? This is slow. default = False'
)


ARGS = CLI.parse_args()

A = ARGS.filename1
B = ARGS.filename2
C = ARGS.sigmathreshold
D = ARGS.plot

#disable printing pandas warnings
pd.options.mode.chained_assignment = None

if (D):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        D = False
        print('Matplotlib not found... not plotting')

## Intro messages
print(' ')
print('Welcome to Flipper!')
time.sleep(1)
print(' ')
print('Brought to you by the Fischer Lab at St. Jude Children\'s Research Hospital')
print('Please cite: ')
print(' ')
time.sleep(1)
print('This program finds peaks in Ringer plots and integrates them.')
print('Flipper detects gain/losses in side chain conformers between ')
print('two datasets and cases where the predominate peak changes ')
print('are FLIPS!')
time.sleep(1)
print(' ')
print(' ')
print('Type -h for info about options')
print(' ')
time.sleep(5)
print('Let\'s get started')
print(' ')
print('Loading the first dataset...')
print('Finding peaks....')
if (D):
    print('Plotting...')

## Part 1: load data, find and integrate peaks, create and export a single dataframe
chis = ['chi1','chi2','chi3','chi4']
df1 = pd.DataFrame([])
for i in chis:
    tmp = peak_find(A,i,C,D)
    df1 = pd.concat([df1,tmp])
print('Peak info from file1 saved to: '+'peak_finder_'+A)
print(' ')
df1.to_csv('peak_finder_'+A,header=True,index=False)
df2 = pd.DataFrame([])
print('Loading the second dataset...')
print('Finding peaks....')
if (D):
    print('Plotting...')
for i in chis:
    tmp = peak_find(B,i,C,D)
    df2 = pd.concat([df2,tmp])
print('Peak info from file2 saved to: '+'peak_finder_'+B)
df2.to_csv('peak_finder_'+B,header=True,index=False)

## Part 2: identify flips
print(' ')
print('Identifying FLIPS....')
print(' ')
for i in chis:
    flips = find_flips(df1,df2,i)
print(' ')

## Part 3: find gain/losses in # of conformers
print('Counting GAIN/LOSS of conformers....')
print('Peaks in each residue in file1 minus file2...')
print('Only printing residues with changes...')
time.sleep(5)
gainloss = pd.DataFrame([])
for i in chis:
    tmp = find_gainloss(df1,df2,i)
    gainloss = pd.concat([gainloss,tmp])
gainloss = gainloss.sort_values(by='res',ascending=True)
gainloss.to_csv(A[:-4]+'_'+B[:-4]+'_gain_loss_peaks.csv',header=True,index=False)
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(gainloss[gainloss['peak_gain_loss']!=0])
print(' ')
print('Peak gain/loss saved to: '+A[:-4]+'_'+B[:-4]+'_gain_loss_peaks.csv')
print(' ')
print('Finished!')
print(' ')


