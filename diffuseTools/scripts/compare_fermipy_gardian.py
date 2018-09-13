#!/usr/bin/env python

import sys

import argparse
import yaml
import numpy as np

from fermipy import utils
from diffuseTools.compare_tools import make_diff_maps, compare_counts


def main():
    usage = "convertGaridanResults.py [options]" 
    description = "Convert Gardian text file results to YAML format used by fermipy"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input Prefix')
    parser.add_argument('-g', '--gardian', type=str, default=None,
                        help='Garidan Prefix ')
    parser.add_argument('-d', '--conv_dict', type=str, default=None,
                        help='Conversion dict ')

    args = parser.parse_args(sys.argv[1:])

    conv_dict = utils.load_yaml(args.conv_dict)

    model_maps_st = []
    model_maps_garg = []

    comp_map = conv_dict['COMP_MAP']
    zmax_map = conv_dict['ZMAX_MAP']
    src_dict = conv_dict['SRC_NAME_DICT']

    do_diff_maps = True
    do_counts = True
    do_map = False

    if do_diff_maps:
        make_diff_maps(comp_map, args.input)

    if do_counts:
        compare_counts(comp_map, args.input, args.gardian, conv_dict)


if __name__ == '__main__':
    main()
