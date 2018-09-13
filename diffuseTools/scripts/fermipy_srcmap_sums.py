#!/usr/bin/env python

import sys

import argparse
import numpy as np
import fermipy.utils as utils
from astropy.io import fits
from fermipy import fits_utils

from diffuseTools.fermipy_tools import get_sums

def main():
    usage = "srcmap_sums.py [options]" 
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
    conv_dict = utils.load_yaml(args.conv_dict)
    comps = conv_dict['COMPONENTS']
    src_dict = conv_dict['SRC_NAME_DICT']

    for c in comps:
        fp = "%s_%s_GAL_V2.fits"%(args.input, c)
        fin = fits.open(fp)
        ebins = fits_utils.find_and_read_ebins(fin)
        for hdu in fin[4:]:
            if hdu.name in ['__weights__']:
                continue
            if hdu.name[0:5] == 'FL8Y ':
                continue
            sums = get_sums(hdu, ebins)
            src_name = src_dict[hdu.name]
            if not out_dict.has_key(src_name):
                out_dict[src_name] = []            
            out_dict[src_name].append(sums)
    
    utils.write_yaml(out_dict, args.output)

if __name__ == '__main__':
    main()
