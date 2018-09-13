#!/usr/bin/env python

import sys

import argparse
import yaml
import numpy as np
from fermipy.utils import load_yaml, write_yaml
from diffuseTools.gardian_tools import make_npred_dict


def main():
    usage = "npred_txt2yaml.py [options]" 
    description = "Convert Gardian text files to yaml dictionary"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output ')
    parser.add_argument('-d', '--conv_dict', type=str, default=None,
                        help='Output ')

    args = parser.parse_args(sys.argv[1:])
    
    conv_dict = load_yaml(args.conv_dict)
    npred_dict = make_npred_dict(args.input, conv_dict)

    write_yaml(npred_dict, args.output)

if __name__ == '__main__':
    main()
