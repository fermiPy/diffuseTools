#!/usr/bin/env python

import sys

import argparse
import yaml
import numpy as np
from fermipy.utils import load_yaml, write_yaml
from diffuseTools.gardian_tools import make_sum_dict

def main():
    usage = "npred_mapSums.py [options]" 
    description = "Convert Gardian text files to yaml dictionary"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output ')
    parser.add_argument('-d', '--conv_dict', type=str, default=None,
                        help='Conversion dict ')
    args = parser.parse_args(sys.argv[1:])

    conv_dict = load_yaml(args.conv_dict)
    comps = conv_dict['COMPONENTS']
    src_dict = conv_dict['SRC_NAME_DICT']
    sel_dep_comps = conv_dict['SEL_DEP_COMPS']
    
    sum_dict = make_sum_dict(args.input, comps, src_dict, sel_dep_comps)
    write_yaml(sum_dict, args.output)


if __name__ == '__main__':
    main()
