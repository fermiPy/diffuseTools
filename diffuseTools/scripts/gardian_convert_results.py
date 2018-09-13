#!/usr/bin/env python

import sys

import argparse
import yaml
import numpy as np

from diffuseTools.gardian_tools import read_lines, get_status, get_var_dict, get_src_dict, translate_src_dict,\
    adjust_xco, adjust_fl8y, pretty_print_input, pretty_print_result


def main():

    usage = "convertGardianResults.py [options]" 
    description = "Convert Gardian text file results to YAML format used by fermipy"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output ')
    parser.add_argument('-c', '--catalog', type=str, default=None,
                        help='Catalog ')
    parser.add_argument('-d', '--conv_dict', type=str, default=None,
                        help='Conversion dictionary ')
    parser.add_argument('--pars', action='store_true', 
                        default=False, help='Write output in form needed to set parameters')
    

    args = parser.parse_args(sys.argv[1:])
    lines = read_lines(args.input)       
    success, mll = get_status(lines)
    var_dict = get_var_dict(lines)
    src_dict = get_src_dict(var_dict)

    if args.conv_dict is not None:
        conv_dict = yaml.load(open(args.conv_dict))
    else:
        print ("You must specify a conversion dictionary with -d option")
        sys.exit(1)

    src_dict_st = translate_src_dict(src_dict, conv_dict)
    adjust_xco(src_dict_st)

    if args.catalog is not None:
        cat_dict = yaml.load(open(args.catalog))
    else:
        cat_dict = None

    src_dict_ad = adjust_fl8y(src_dict_st, cat_dict)

    if args.pars:
        pretty_print_input(src_dict_ad, args.output)
    else:
        pretty_print_result(src_dict_ad, args.output)


if __name__ == '__main__':
    main()
