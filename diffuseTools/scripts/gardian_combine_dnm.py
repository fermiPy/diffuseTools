#!/usr/bin/env python

import sys
import os
import argparse

import numpy as np
import fermipy.utils as utils


from diffuseTools.gardian_tools import combine_dnm_maps


def main():
    usage = "combine_dnm.py [options]" 
    description = "add the pos and neg dnm components"

    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file prefix')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output ')
    parser.add_argument('-f', '--fact', default=False,
                        action='store_true', help='Use factors')    
    
    args = parser.parse_args(sys.argv[1:])

    out_dict = utils.load_yaml(args.input)
    
    combine_dnm_maps(out_dict, args.fact)

    utils.write_yaml(out_dict, args.output)


if __name__ == '__main__':
    main()
