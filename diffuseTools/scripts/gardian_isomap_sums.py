#!/usr/bin/env python

import os
import sys

import argparse
import yaml
import numpy as np
from fermipy.utils import write_yaml

from diffuseTools.gardian_tools import get_expmap_sums, get_int_egb

COMPONENTS = ['E0_PSF3',
              'E1_PSF23',
              'E2_PSF123',
              'E3_PSF0123']

def main():
    usage = "npred_mapSums.py [options]" 
    description = "Convert Gardian text files to yaml dictionary"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file')
    parser.add_argument('-e', '--exposure', type=str, default=None,
                        help='Exposure file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output ')
    args = parser.parse_args(sys.argv[1:])

    erange_0_edges = np.logspace(np.log10(30),2,3)
    erange_1_edges = np.logspace(2,np.log10(300),3)
    erange_2_edges = np.logspace(np.log10(300),3,4)
    erange_3_edges = np.logspace(3,6,10)   
        
    erange_0 = np.sqrt(erange_0_edges[0:-1]*erange_0_edges[1:])
    erange_1 = np.sqrt(erange_1_edges[0:-1]*erange_1_edges[1:])
    erange_2 = np.sqrt(erange_2_edges[0:-1]*erange_2_edges[1:])
    erange_3 = np.sqrt(erange_3_edges[0:-1]*erange_3_edges[1:])

    energies = [erange_0, erange_1, erange_2, erange_3]
    ewidths = [erange_0_edges[1:] - erange_0_edges[0:-1],
               erange_1_edges[1:] - erange_1_edges[0:-1],
               erange_2_edges[1:] - erange_2_edges[0:-1],
               erange_3_edges[1:] - erange_3_edges[0:-1]]


    l = []
    factor = 12*32*32
    for evals, ew, comp in zip(energies, ewidths, COMPONENTS):
        exppath = os.path.join(args.exposure,"%s_expcube_moreplanes_ring_P8R3.fits"%(comp))
        egbpath = args.input

        exp_vals = get_expmap_sums(exppath, evals)
        egb_vals = get_int_egb(egbpath, evals)
        
        o_vals = exp_vals*egb_vals*ew/factor
        l.append(o_vals)

    o = dict(iso_map=l)
    write_yaml(o, args.output)


if __name__ == '__main__':
    main()
