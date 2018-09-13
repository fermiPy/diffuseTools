# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools to interact with matplotlib
"""
from __future__ import absolute_import, division, print_function

import sys

import argparse
import yaml
import numpy as np
import healpy as hp

import pickle
import fermipy.gtutils as gtutils
import fermipy.utils as utils
import pyLikelihood as pyLike
import matplotlib.pyplot as plt
from matplotlib.colors import ColorConverter
from astropy.io import fits

FIGURE_LAYOUT = {1:(1,1),
                 2:(1,2),
                 3:(2,2),
                 4:(2,2),
                 5:(2,3),
                 6:(2,3),
                 7:(3,3),
                 8:(3,3),
                 9:(3,3),
                 10:(3,4),
                 11:(3,4),
                 12:(3,4)}

CC = ColorConverter()
GREY = CC.to_rgba('grey')
BLUE = CC.to_rgba('blue')

def plot_map_sums(map_dict, srclist, compare_map_dict=None):
    
    (nrows, ncols) = FIGURE_LAYOUT[len(srclist)]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14,8))
    try: 
        energies = map_dict['Energies']
    except KeyError as msg:
        print (map_dict.keys())
        raise KeyError(msg)
    if compare_map_dict is not None:
        energies_comp = compare_map_dict['Energies']
        e_widths = compare_map_dict['E_widths']
    for i,k in enumerate(srclist):
        irow = int(i/ncols)
        icol = int(i - ncols*(irow))
        if nrows > 1 and ncols > 1:
            ax = axs[irow, icol]
        elif nrows > 1:
            ax = axs[irow]
        elif ncols > 1:
            try:
                ax = axs[icol]
            except IndexError as msg:
                print (icol)
                raise IndexError(msg)
        else:
            ax = axs
        ax.set_title(k)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('E [MeV]')
        ax.set_ylabel('Spectrum [a.u.]')   
        ax.set_xlim(10, 1e6)
        try:
            v_list = map_dict[k]
            for v,e in zip(v_list, energies):
                try:
                    ax.plot(e, v, color=GREY)
                except ValueError, msg:
                    print (msg)
        except KeyError:
            pass
        if compare_map_dict is not None:
            try:
                vc_list = compare_map_dict[k]
                for vc,e,w in zip(vc_list, energies_comp, e_widths):
                    
                    try:
                        ax.plot(e, vc/np.array(w), color=BLUE)
                    except ValueError, msg:
                        print (msg)
            except KeyError, msg:
                print (msg)
    return fig


def plot_map_diff(map_dict, srclist, compare_map_dict):
    
    (nrows, ncols) = FIGURE_LAYOUT[len(srclist)]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14,8))
    energies = compare_map_dict['Energies']
    e_widths = compare_map_dict['E_widths']
    for i,k in enumerate(srclist):
        irow = int(i/ncols)
        icol = int(i - ncols*(irow))
        if nrows > 1 and ncols > 1:
            ax = axs[irow, icol]
        elif nrows > 1:
            ax = axs[irow]
        elif ncols > 1:
            ax = axs[icol]
        else:
            ax = axs
        ax.set_title(k)
        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.set_xlabel('E [MeV]')
        ax.set_ylabel('Spectrum [a.u.]')   
        ax.set_xlim(10, 1e6)
        try:
            v_list = map_dict[k]
            vc_list = compare_map_dict[k]
            for v,vc,w,e in zip(v_list, vc_list, e_widths, energies):
                va = np.array(v)
                vm = np.sqrt(va[1:]*va[0:-1])
                frac_diff = (vm - (np.array(vc)/np.array(w)))/vm
                try:
                    ax.plot(e, frac_diff, color='b')
                except ValueError, msg:
                    print (msg)
        except KeyError:
            pass
    return fig



def sum_npreds(npred_dict):
    rsum = None
    for k,v in npred_dict.items():
        if rsum is None:
            rsum = np.array(v)
        else:
            rsum += np.array(v)
    npred_dict['total'] = rsum


def plot_npreds(npred_dict, srclist, compare_npred_dict=None):
    
    (nrows, ncols) = FIGURE_LAYOUT[len(srclist)]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14,8))
    energies = npred_dict['Energies']
    if compare_npred_dict is not None:
        energies_comp = compare_npred_dict['Energies']
    
    for i,k in enumerate(srclist):
        irow = int(i/ncols)
        icol = int(i - ncols*(irow))
        if nrows > 1 and ncols > 1:
            ax = axs[irow, icol]
        elif nrows > 1:
            ax = axs[irow]
        elif ncols > 1:
            ax = axs[icol]
        else:
            ax = axs
        ax.set_title(k)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('E [MeV]')
        ax.set_ylabel('Spectrum [a.u.]')   
        ax.set_xlim(10, 1e6)
        try:
            v = npred_dict[k]
            try:
                ax.plot(energies, v, color=GREY)
            except ValueError:
                pass
        except KeyError:
            pass
        if compare_npred_dict is not None:
            try:
                vc_list = compare_npred_dict[k]
                for vc,e in zip(vc_list, energies_comp):
                    try:
                        ax.plot(e, vc, color=BLUE)
                    except ValueError:
                        pass
            except KeyError:
                pass        
    return fig


def plot_npred_diff(npred_dict, srclist, compare_npred_dict=None):
    
    (nrows, ncols) = FIGURE_LAYOUT[len(srclist)]
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(14,8))
    energies = np.hstack(compare_npred_dict['Energies'])
    
    for i,k in enumerate(srclist):
        irow = int(i/ncols)
        icol = int(i - ncols*(irow))
        if nrows > 1 and ncols > 1:
            ax = axs[irow, icol]
        elif nrows > 1:
            ax = axs[irow]
        elif ncols > 1:
            ax = axs[icol]
        else:
            ax = axs
        ax.set_title(k)
        ax.set_xscale('log')
        ax.set_yscale('linear')
        ax.set_xlabel('E [MeV]')
        ax.set_ylabel('Spectrum [a.u.]')   
        ax.set_xlim(10, 1e6)
        try:
            v = np.array(npred_dict[k])
            vc_list = compare_npred_dict[k]
            vc = np.hstack(vc_list)
            #total = npred_dict['total']
            total = v + vc
            frac_diff = (v - vc)/total
            mask = vc > 0
            try:
                ax.plot(energies[mask], frac_diff[mask], color='b')
            except ValueError:
                pass
        except KeyError:
            pass      
    return fig



def make_grid():
    l2 = [-90, 60, 60, -90, -90]
    b2 = [-70, -70, 80, 80, -70]
    larr = []
    barr = []
    numpts = 100
    for k in range(len(l2)-1):
        slopex = (l2[k+1] - l2[k])/float(numpts)
        slopey = (b2[k+1] - b2[k])/float(numpts)
        larr.extend(l2[k] + np.array(xrange(numpts))*slopex)
        barr.extend(b2[k] + np.array(xrange(numpts))*slopey)
        
    theta = np.pi/2 - np.array(barr)*np.pi/180
    phi = np.array(larr)*np.pi/180
    return (theta, phi)



def make_regions():
    regions = [ ( (0,  360), (-20,20) ),
                ( (300, 60), (-10,10) ),
                ( (0,  360), (60,90) ),
                ( (0,  360), (-90,-60) ),
                ( (0,  360), (-40,-20) ),
                ( (0,  360), (20,40) ),
                ( (98, 148), (40,60) ),
                ( (278, 328), (-60,-40) ),
                ( (98, 148), (-60,-40) ),
                ( (278, 328), (40,60) ),
                ( (0,  360), (-90,90) ),
                ]
    
    regionNames = [ "GalPlane",
                    "GalRidge",
                    "NorthPole",
                    "SouthPole",
                    "SouthInterm",
                    "NorthInterm",
                    "NorthLimb",
                    "SouthLimb",
                    "SouthAntiLimb",
                    "NorthAntiLimb",
                    "AllSky",]

    pickleFile = '/gpfs/slac/kipac/fs1/u/dmcat/data/flight/diffuse_fitting/analysis/scripts/regionPixels.pickle'
    infile = open(pickleFile, 'rb')
    regionPixels = pickle.load(infile)
    infile.close()
    regionPixels.append(None)
    
    return dict(regions=regions,
                regionNames=regionNames,
                regionPixels=regionPixels)
    

def figure_model(model, erange):    
    
    temp = np.sort(model)
    minv = temp[int(0.005*len(model))]
    maxv = temp[int(0.995*len(model))]

    hp.mollview(model,
                norm='log',
                title='Model %.2f-%.2f GeV'%(erange[0], erange[1]),
                sub=[1,3,1],
                min=minv, max=maxv,
                cmap='jet')
    


def figure_frac_resid(counts, model, erange):
    
    frac_resid = (counts - model)/model
    temp = np.sort(frac_resid)
    minv = temp[int(0.005*len(frac_resid))]
    maxv = temp[int(0.995*len(frac_resid))]
    hp.mollview(frac_resid,
                title='Frac. Residual %.2f-%.2f GeV'%(erange[0], erange[1]),
                sub=[1,3,2],
                min=minv, max=maxv,
                cmap='jet')

    


def figure_signif(counts, model, erange):

    GRID = make_grid()

    signif = (counts - model)/np.sqrt(model)
    minv = -5.
    maxv = 5.
    hp.mollview(signif,
                title='Residual Signif. %.2f-%.2f GeV'%(erange[0], erange[1]),
                sub=[1,3,3],
                min=minv, max=maxv,
                cmap='jet')

    hp.projplot(GRID[0], GRID[1], 'b,')
    


def figure_region_spectra(regionNames, regionSpectra, ecents, ewidths):
    
    row = 3
    col = 4
    fig, axs = plt.subplots(nrows=row, ncols=col, sharex=True, figsize=(14,8))
    
    for i, regName in enumerate(regionNames):

        
        counts = regionSpectra[regName]['counts'] * ecents
        model = regionSpectra[regName]['model'] * ecents
        ax = axs[int(i/4), i - 4*int(i/4)]
        ax.set_xscale('log')
        ax.set_yscale('log')

        ymin = min( counts.min(), model.min() ) * 0.1
        ymax = min( counts.max(), model.max() ) * 10.

        ax.errorbar(ecents, np.maximum(1, counts), xerr=ewidths/2., label='counts')
        ax.errorbar(ecents, np.maximum(1, model), xerr=ewidths/2., label='model')
        ax.set_ylim([ymin, ymax])

        axs[2,int(i/row)].set_xlabel("E [MeV]")
        axs[int(i/col),0].set_ylabel("Counts*E per bin")
        
        ax.set_title(regName,fontsize=12)
        legend = ax.legend(loc='upper right',prop={'size':6})
        for label in legend.get_texts():
            label.set_fontsize('8')

    return fig


def figure_region_resid(regionNames, regionSpectra, ecents, ewidths):
    
    row = 3
    col = 4
    fig, axs = plt.subplots(nrows=row, ncols=col, sharex=True, figsize=(14,8))
    
    for i, regName in enumerate(regionNames):
        counts = regionSpectra[regName]['counts']
        model = regionSpectra[regName]['model']

        frac_resid = (counts - model) / model
        ax = axs[int(i/4), i - 4*int(i/4)]
        ymin = min(frac_resid) - 0.1
        ymax = max(frac_resid) + 0.1        
        
        ax.semilogx(ecents, frac_resid)
        ax.set_ylim([ymin, ymax])
        axs[2,int(i/row)].set_xlabel("E [MeV]")
        axs[int(i/col),0].set_ylabel("(count - model)/model",fontsize=10)
        ax.set_title(regName,fontsize=12)
        ax.axhline(0, color='black', ls='--')

    return fig


