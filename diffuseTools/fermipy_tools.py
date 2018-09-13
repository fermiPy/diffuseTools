# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools to interact with fermipy
"""
from __future__ import absolute_import, division, print_function

import sys

import numpy as np
from fermipy.gtanalysis import GTAnalysis
from fermipy.skymap import HpxMap


def build_srcdict(gta, prop):
    """ Build a dictionary that maps src to a property """
    o = {}
    for s in gta.roi.sources:
        o[s.name] = s[prop]
    return o

def get_src_names(gta):
    o = []
    for s in gta.roi.sources:
        o += [s.name]
    return sorted(o)    

def set_wts_get_npred_wt(gta, maskname):
    gta.set_weights_map(maskname, update_roi=False)
    sys.stdout.write("Updating %i sources: "%(len(gta.like.sourceNames())))
    sys.stdout.flush()
    for i,name in enumerate(gta.like.sourceNames()):
        if i%20==0:
            sys.stdout.write('x')
        else:
            sys.stdout.write('.')
        sys.stdout.flush()
        gta._init_source(name)
    gta._update_roi()
    sys.stdout.write('!\n')
    return build_srcdict(gta, 'npred_wt')




def make_npred_dict(src_dict):
    o = {}
    for k,v in sorted(src_dict.items()):
        o[k] = v['model_counts_wt']
    return o

def get_energies(roi_dict):
    ebins = roi_dict['energies']
    evals = np.sqrt(ebins[0:-1]*ebins[1:])    
    return evals


def get_sums(hdu, ebins):
    try:
        themap = HpxMap.create_from_hdu(hdu, ebins)
    except:
        print ("Failed to read %s"%(hdu.name))
        return np.zeros((ebins.shape))
    return themap.counts.sum(1)
