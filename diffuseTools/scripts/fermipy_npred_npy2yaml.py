#!/usr/bin/env python

import sys

import argparse
import yaml
import numpy as np
from fermipy.utils import write_yaml, load_npy
from diffuseTools.fermipy_tools import make_npred_dict, get_energies


def main():
    usage = "npred_txt2yaml.py [options]" 
    description = "Convert Gardian text files to yaml dictionary"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output ')
    args = parser.parse_args(sys.argv[1:])

    src_dict = load_npy(args.input)
    npred_dict = make_npred_dict(src_dict['sources'])
    npred_dict['Energies'] = get_energies(src_dict['roi'])

    write_yaml(npred_dict, args.output)

if __name__ == '__main__':
    main()
