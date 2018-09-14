# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools to compare fermipy and GARDIAN results
"""
from __future__ import absolute_import, division, print_function

import sys

import argparse
import yaml
import numpy as np

from gammapy.maps import Map, HpxNDMap
from fermipy.utils import load_npy

from diffuseTools.gardian_tools import read_lines, read_comp_names, read_values, read_gardian_counts


def compare_counts_dicts(st_counts, gard_counts, ostream=sys.stdout):
    for k, v in sorted(st_counts.items()):
        try:
            vv = gard_counts.pop(k)
        except KeyError:
            ostream.write("%-25s : %.4e : ST ONLY \n"%(k, v.sum()))
            continue
        diff = v.sum() - vv.sum()
        fdiff = diff / ( v.sum() + vv.sum() )
        ostream.write("%-25s : %.4e : %.4e : %.4e : %.4e \n"%(k, v.sum(), vv.sum(), diff, fdiff))
    for k, v in sorted(gard_counts.items()):
        ostream.write("%-25s :  GA only   : %.4e\n"%(k, v.sum()))


def make_diff_maps(comp_map, input_pref):

    for k,v in sorted(comp_map.items()):
        fname_st = "mcube_%s_%s.fits"%(input_pref,k)
        #fname_gard = "%s_FinalModels/%s.fits.gz"%(args.gardian, v)
        fname_gard = "../../GardianResults_local/%s.fits.gz"%(v)
        map_st = Map.read(fname_st, 'SKYMAP')
        try:
            map_gard = Map.read(fname_gard, 'SKYMAP2')
        except KeyError as msg:
            print ("Failed to read %s " % fname_gard)
            raise KeyError(msg)
                
        map_gard_cast = map_gard.to_ud_graded(map_st.geom.nside, True)
        diff = map_st.data - map_gard_cast.data
        map_diff = HpxNDMap(map_st.geom, diff)
        map_diff.write("diff_gard_%s.fits"%(k))

        esum = diff.sum(0)
        stsum = map_st.data.sum(0)
        fsum = esum / stsum
        frac_diff = HpxNDMap(map_st.geom.to_image(), fsum)
        frac_diff.write("diff_gard_frac_%s.fits"%(k))


def compare_counts(comp_map, input_pref, gard_pref, src_dict):
    gard_counts = {}
    for k,v in sorted(comp_map.items()):
        fname_gard = "%s_ModelComponentsCounts/%s.txt"%(gard_pref, v)
        counts = read_gardian_counts(fname_gard)
        for kk,vv in sorted(counts.items()):
            if kk in ['Energies', 'E_widths']:
                continue
            src = src_dict[kk]
            if not gard_counts.has_key(src):
                gard_counts[src] = []
                gard_counts[src] += [vv]

    for k,v in sorted(gard_counts.items()):
        gard_counts[k] = np.hstack(v)
    
    st_counts = {}
    dd = load_npy("%s.npy"%(input_pref))['sources']
    for k,v in sorted(dd.items()):
        if k[0:4] == 'FL8Y':
            src = 'FL8Y'
        else:
            src = k
        st_counts[src] = v['model_counts_wt']

    compare_counts_dicts(st_counts, gard_counts)

