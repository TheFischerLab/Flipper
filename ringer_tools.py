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
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def peak_find(file,chi,sigma_threshold,plot,height,prominence,width,distance):
    """ Find peaks in Ringer plots """

    dataframe = pd.read_csv(file,header=None)

    # replaces datapoint with angles from ringer output CSV
    step=5
    angles = np.arange(0,360,step=step)

    output = []

    for i in range(len(dataframe)):
        if dataframe.iloc[i,1] == chi:
            res_n = dataframe.iloc[i,0][6:-1].strip()
            chain = dataframe.iloc[i,0][4]
            res = dataframe.iloc[i,0]
            chi = dataframe.iloc[i,1]
            sigma_raw = dataframe.iloc[i,3:]
            sigma_raw = np.array(sigma_raw)
            dat = np.column_stack((angles,sigma_raw))
            dat = [(x,y) for (x,y) in dat if y >= sigma_threshold]
            sigma_above_threshold = [y for (x,y) in dat]
            angles_cut = [x for (x,y) in dat]
            peaks, properties = find_peaks(sigma_above_threshold,
                                           height=height,
                                           prominence=prominence,
                                           width=width,
                                           distance=distance)
            peak_n = len(peaks)
            area = []
            lefts = []
            rights = []
            for peak in range(peak_n):
                left = properties['left_bases'][peak]
                lefts.append(left)
                right = properties['right_bases'][peak]
                rights.append(right)
                data_to_int = np.column_stack([angles_cut[left:right],
                                               sigma_above_threshold[left:right]])
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
                    for areai in area:
                        areas_norm.append(areai/sum_area)
                else:
                    sum_area = sum(area)
                    for areai in area:
                        areas_norm.append(areai/sum_area)
            if len(area) == 1:
                areas_norm = 1.0
            if len(area) == 0:
                areas_norm = 0
            if plot:
                plt.figure(figsize=(5,5))
                plt.plot(angles,sigma_raw,color='black')
                plt.plot(angles_cut,sigma_above_threshold,'red',linestyle='dotted')
                plt.plot([angles_cut[x] for x in peaks],
                         [sigma_threshold]*len(peaks), "x")
                plt.axhline(y=sigma_threshold, linestyle="--", color="gray")
                plt.legend(['raw data',
                            'data used for peak \n finding and integration',
                            'peak','sigma cutoff'],
                            bbox_to_anchor=(1.05, 1),
                            loc='upper left')
                plt.xlabel(chi+'  angle (˚)')
                plt.ylabel('sigma')
                plt.title(res)
                if len(peaks)>0:
                    plt.ylim(np.min(sigma_raw),np.max(sigma_raw)+sigma_threshold)
                plt.tight_layout()
                plt.savefig(str(res_n)+'_'+chain+'_'+file+'_'+chi+'_peaks.png',dpi=100)
                plt.close()
            peaks = [angles_cut[x] for x in peaks]
            output.append((res_n,res,chi,peaks,peak_n,areas_norm,chain))
    output = pd.DataFrame(output)
    output.columns = ['res_n','res','chi','peak_angles','peaks_n','areas_norm','chain']
    return output


def find_flips(peaks1,peaks2,chi):
    """Find changes (flips) in major/minor conformation of a residue between
       two datafiles based on peak integration """
    dataframe_f1 = peaks1
    dataframe_f2 = peaks2

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
    dataframe_f1 = dataframe_f1[dataframe_f1['peaks_n']==2]
    dataframe_f2 = dataframe_f2[dataframe_f2['peaks_n']==2]

    # only consider one chi angle
    dataframe_f1 = dataframe_f1[dataframe_f1['chi']==chi]
    dataframe_f2 = dataframe_f2[dataframe_f2['chi']==chi]

    # only consider residues in allowable chi list
    dataframe_f1['res'] = dataframe_f1['res'].str[:3]
    dataframe_f1 = dataframe_f1[dataframe_f1['res'].isin(good_residues)]
    dataframe_f2['res'] = dataframe_f2['res'].str[:3]
    dataframe_f2 = dataframe_f2[dataframe_f2['res'].isin(good_residues)]

    # format
    dataframe_f1['areas_norm'] = dataframe_f1['areas_norm'].astype(str).str[1:-1]
    dataframe_f1['areas_norm'] = [x.split(',') for x in dataframe_f1['areas_norm']]
    dataframe_f1['areas_norm'] = [(float(x),float(y)) for x,y in dataframe_f1['areas_norm']]
    dataframe_f2['areas_norm'] = dataframe_f2['areas_norm'].astype(str).str[1:-1]
    dataframe_f2['areas_norm'] = [x.split(',') for x in dataframe_f2['areas_norm']]
    dataframe_f2['areas_norm'] = [(float(x),float(y)) for x,y in dataframe_f2['areas_norm']]

    #format
    dataframe_f1['res_n'] = dataframe_f1['res_n'].astype(int)
    dataframe_f2['res_n'] = dataframe_f2['res_n'].astype(int)

    flips_pair = []

    for j in range(len(dataframe_f1)):
        res_n = int(dataframe_f1.iloc[j,0])
        if res_n in dataframe_f2['res_n'].values:
            f1_p1 = dataframe_f1.iloc[j,5][0]
            f2_p1 = dataframe_f2[dataframe_f2['res_n']==res_n]['areas_norm'].values[0][0]
            if f1_p1 > 0.5 > f2_p1:
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
            if f1_p1 < 0.5 < f2_p1:
                print('FLIP! at residue',res_n,chi)
                flips_pair.append((res_n,chi))
    flips = pd.DataFrame(flips_pair,columns=['res_n',chi])
    return flips

def find_gainloss(peaks1,peaks2,chi):
    """ Calculate the changes in number of peaks (conformations) between two datafiles """
    dataframe_f1 = peaks1
    dataframe_f2 = peaks2

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
    dataframe_f1 = dataframe_f1[dataframe_f1['chi']==chi]
    dataframe_f2 = dataframe_f2[dataframe_f2['chi']==chi]

    # only consider residues in allowable chi list
    dataframe_f1['res'] = dataframe_f1['res'].str[:3]
    dataframe_f1 = dataframe_f1[dataframe_f1['res'].isin(good_residues)]
    dataframe_f2['res'] = dataframe_f2['res'].str[:3]
    dataframe_f2 = dataframe_f2[dataframe_f2['res'].isin(good_residues)]

    dataframe_f1 = dataframe_f1[['res_n','peaks_n','res']]
    dataframe_f1.columns = ['res_n','peak_cryo','res']
    dataframe_f2 = dataframe_f2[['res_n','peaks_n','res']]
    dataframe_f2.columns = ['res_n','peak_rt','res']

    dataframe = pd.merge(dataframe_f1,dataframe_f2,on='res_n',how='left')
    dataframe = dataframe.fillna(0)

    differences = dataframe.peak_cryo.values - dataframe.peak_rt.values
    residues = dataframe['res_n'].values[:]

    dataframe = np.column_stack((residues,differences))

    dataframe = pd.DataFrame(dataframe,columns = ['res','peak_gain_loss'])
    dataframe['res'] = dataframe['res'].astype(int)
    dataframe['chi'] = chi

    return dataframe
