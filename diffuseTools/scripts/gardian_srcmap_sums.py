#!/usr/bin/env python

import sys
import os

import argparse
import numpy as np
import fermipy.utils as utils
from astropy.io import fits
from fermipy.skymap import read_map_from_fits
from fermipy import fits_utils


def main():
    usage = "fermitool-gardian_srcmap_sums.py [options]" 
    description = "sum source maps files over all pixels"

    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file prefix')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output ')
    parser.add_argument('-d', '--conv_dict', type=str, default=None,
                        help='Conversion dict ')
    args = parser.parse_args(sys.argv[1:])

    out_dict = {}    

    conv_dict = load_yaml(args.conv_dict)
    comps = conv_dict['COMPONENTS']
    src_dict = conv_dict['SRC_NAME_DICT']
    sel_dep_comps = conv_dict['SEL_DEP_COMPS']
    cache_dirs = conv_dict['CACHE_DIR_DICT']

    sys.stdout.write("Reading %i components\n"%(len(SRC_NAME_DICT)))
    
    for k,v in sorted(src_dict.items()):
        sys.stdout.write('.')
        sys.stdout.flush()
        cache_dir = cache_dirs[k]
        for c in comps:
            if k in sel_dep_comps:
                fp = os.path.join(args.input,cache_dir, "%s_%s%s_cache.fits.gz"%(c, k, c[0:2]))
            else:
                fp = os.path.join(args.input,cache_dir, "%s_%s_cache.fits.gz"%(c, k))

            map_in = read_map_from_fits(fp)        
            sums = map_in.counts.sum(1)            
            if not out_dict.has_key(v):
                out_dict[v] = []            
            out_dict[v].append(sums)

    sys.stdout.write('!\n')
    
    utils.write_yaml(out_dict, args.output)


if __name__ == '__main__':
    main()
