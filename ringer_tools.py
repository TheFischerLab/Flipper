#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Tools for analyzing Ringer measurements. 
#    Part of the Flipper program.
#    Authors: Tim Stachowski & Marcus Fischer
#    Email: tim.stachowski@stjude.org
#    Copyright 2022 St. Jude Children's Research Hospital, Memphis, TN.
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

import numpy as np
import pandas as pd
from itertools import product
from scipy.signal import find_peaks
import argparse
import time

def peak_find(file,chi,sigma_threshold,plot,height,prominence,width,distance):
    df = pd.read_csv(file,header=None)

    # replaces datapoint with angles from ringer output CSV
    step=5
    angles = np.arange(0,360,step=step)
    columns = ['res','chi','rand']+list(angles)

    chi = chi
    sigma_threshold = sigma_threshold
    plot = plot
    output = []

    for i in range(len(df)):
        if df.iloc[i,1] == chi:
            res_n = df.iloc[i,0][6:-1].strip()
            chain = df.iloc[i,0][4]
            res = df.iloc[i,0]
            chi = df.iloc[i,1]
            sigma_raw = df.iloc[i,3:]
            sigma_raw = np.array(sigma_raw)
            dat = np.column_stack((angles,sigma_raw))
            dat = [(x,y) for (x,y) in dat if y >= sigma_threshold]
            sigma_above_threshold = [y for (x,y) in dat]
            angles_cut = [x for (x,y) in dat]
            peaks, properties = find_peaks(sigma_above_threshold, height=sigma_threshold,prominence=prominence,width=width,distance=distance)
            peak_n = len(peaks)
            area = []
            lefts = []
            rights = []
            for i in range(peak_n):
                left = properties['left_bases'][i]
                lefts.append(left)
                right = properties['right_bases'][i]
                rights.append(right)
                data_to_int = np.column_stack([angles_cut[left:right],sigma_above_threshold[left:right]])
                angles_int = [j for (j,y) in data_to_int]
                sigma_int = [y for (j,y) in data_to_int]
                area.append(np.trapz(sigma_int,angles_int))
            areas_norm = []
            if len(area) > 1:
                if (lefts[0] == lefts[-1]) or (rights[0] == rights[-1]):
                    if area[0] > area[-1]:
                        area[0] = area[0]-area[-1]
                        area[0] = np.absolute(area[0])
                    else:
                        area[-1] = area[-1] - area[0]
                        area[-1] = np.absolute(area [-1])
                    sum_area = sum(area)
                    for i in area:
                        areas_norm.append(i/sum_area)
                else:
                    sum_area = sum(area)
                    for i in area:
                        areas_norm.append(i/sum_area)
            if len(area) == 1:
                areas_norm = 1.0
            if len(area) == 0:
                areas_norm = 0
            if (plot):
                plt.figure(figsize=(5,5))
                plt.plot(angles,sigma_raw,color='black')
                plt.plot(angles_cut,sigma_threshold,'red',linestyle='dotted')
                plt.plot([angles_cut[x] for x in peaks], [threshold]*len(peaks), "x")
                plt.axhline(y=threshold, linestyle="--", color="gray")
                plt.legend(['raw data','data used for peak \n finding and integration','peak','sigma cutoff'],bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.xlabel(chi+'  angle (˚)')
                plt.ylabel('sigma')
                plt.title(res)
                if len(peaks)>0:
                    plt.ylim(np.min(sigma_raw),np.max(sigma_raw)+threshold)
                plt.tight_layout()
                plt.savefig(str(res_n)+'_'+chain+'_'+file+'_'+chi+'_peaks.png',dpi=100)
                plt.close()
            peaks = [angles_cut[x] for x in peaks]
            output.append((res_n,res,chi,peaks,peak_n,areas_norm,chain))
    output = pd.DataFrame(output)
    output.columns = ['res_n','res','chi','peak_angles','peaks_n','areas_norm','chain']
    return output

def find_flips(peaks1,peaks2,chi):
    
    df_cryo = peaks1
    df_rt = peaks2

    if chi == 'chi1':
        good_residues = "Ser, Gln, Asn, Glu, Asp, Arg, Lys, Met, Cys, Leu"
    if chi == 'chi2':
        good_residues = "Gln, Glu, Arg, Lys, Met, Ile"
    if chi == 'chi3':
        good_residues = "Lys, Arg, Met"
    if chi == 'chi4':
        good_residues = "Lys, Arg"
    good_residues = good_residues.split(', ')
    good_residues = [x.upper() for x in good_residues]
    
    # only consider ringer plots with 2 peaks
    df_cryo = df_cryo[df_cryo['peaks_n']==2]
    df_rt = df_rt[df_rt['peaks_n']==2]
    
    # only consider one chi angle
    df_cryo = df_cryo[df_cryo['chi']==chi]
    df_rt = df_rt[df_rt['chi']==chi]

    # only consider residues in allowable chi list 
    df_cryo['res'] = df_cryo['res'].str[:3]
    df_cryo = df_cryo[df_cryo['res'].isin(good_residues)]
    df_rt['res'] = df_rt['res'].str[:3]
    df_rt = df_rt[df_rt['res'].isin(good_residues)]
    
    # format
    df_cryo['areas_norm'] = df_cryo['areas_norm'].astype(str).str[1:-1]
    df_cryo['areas_norm'] = [x.split(',') for x in df_cryo['areas_norm']]
    df_cryo['areas_norm'] = [(float(x),float(y)) for x,y in df_cryo['areas_norm']]
    df_rt['areas_norm'] = df_rt['areas_norm'].astype(str).str[1:-1]
    df_rt['areas_norm'] = [x.split(',') for x in df_rt['areas_norm']]
    df_rt['areas_norm'] = [(float(x),float(y)) for x,y in df_rt['areas_norm']]
    
    #format
    df_cryo['res_n'] = df_cryo['res_n'].astype(int)
    df_rt['res_n'] = df_rt['res_n'].astype(int)
    
    flips_pair = []
    
    for j in range(len(df_cryo)):
        res_n = int(df_cryo.iloc[j,0])
        if (res_n in df_rt['res_n'].values):
            if (df_cryo.iloc[j,5][0] > 0.5) and (df_rt[df_rt['res_n']==res_n]['areas_norm'].values[0][0] < 0.5):
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
            if (df_cryo.iloc[j,5][0] < 0.5) and (df_rt[df_rt['res_n']==res_n]['areas_norm'].values[0][0] > 0.5):
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
    flips = pd.DataFrame(flips_pair,columns=['res_n',chi])
    return flips

def find_gainloss(peaks1,peaks2,chi):
    
    df_cryo = peaks1
    df_rt = peaks2
    
    if chi == 'chi1':
        good_residues = "Ser, Gln, Asn, Glu, Asp, Arg, Lys, Met, Cys, Leu"
    if chi == 'chi2':
        good_residues = "Gln, Glu, Arg, Lys, Met, Ile"
    if chi == 'chi3':
        good_residues = "Lys, Arg, Met"
    if chi == 'chi4':
        good_residues = "Lys, Arg"
    good_residues = good_residues.split(', ')
    good_residues = [x.upper() for x in good_residues]
    
    # only consider one chi angle
    df_cryo = df_cryo[df_cryo['chi']==chi]
    df_rt = df_rt[df_rt['chi']==chi]

    # only consider residues in allowable chi list 
    df_cryo['res'] = df_cryo['res'].str[:3]
    df_cryo = df_cryo[df_cryo['res'].isin(good_residues)]
    df_rt['res'] = df_rt['res'].str[:3]
    df_rt = df_rt[df_rt['res'].isin(good_residues)]
    
    df_cryo = df_cryo[['res_n','peaks_n','res']]
    df_cryo.columns = ['res_n','peak_cryo','res']
    df_rt = df_rt[['res_n','peaks_n','res']]
    df_rt.columns = ['res_n','peak_rt','res']
    
    df = pd.merge(df_cryo,df_rt,on='res_n',how='left')
    df = df.fillna(0)
    
    differences = df.peak_cryo.values - df.peak_rt.values
    residues = df['res_n'].values[:]
    
    df = np.column_stack((residues,differences))
    
    df = pd.DataFrame(df,columns = ['res','peak_gain_loss'])
    df['res'] = df['res'].astype(int)
    df['chi'] = chi
    
    return df