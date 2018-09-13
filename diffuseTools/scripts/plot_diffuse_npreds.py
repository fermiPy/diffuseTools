#!/usr/bin/env python

import sys

import argparse
import yaml
import numpy as np

import fermipy.utils as utils
from diffuseTools.plot_tools import plot_npreds, plot_npred_diff, sum_npreds

def main():

    usage = "make_spectra.py [options]" 
    description = "Plot spectral components from a results file"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file')    
    parser.add_argument('-c', '--compare', type=str, default=None,
                        help='Input file to compare')    
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output file prefix ')
    parser.add_argument('-p', '--plot', type=str, default='png',
                        help='Plot type')
    parser.add_argument('-d', '--comp_dict', type=str, default=None,
                        help='Component dict')
    args = parser.parse_args(sys.argv[1:])

    if args.input is not None:
        npred_dict = utils.load_yaml(args.input)
    else:
        npred_dict = None

    if args.compare is not None:
        compare_npred_dict = utils.load_yaml(args.compare)
    else:
        compare_npred_dict = None

    comp_dict = utils.load_yaml(args.comp_dict)

    sum_npreds(npred_dict)
    
    figs = []
    for k, v in comp_dict.items():
        fig_map = plot_npreds(npred_dict, v, compare_npred_dict)
        fig_diff = plot_npred_diff(npred_dict, v, compare_npred_dict)
        fig_map.savefig("%s_%s_npred.%s"%(args.output, k, args.plot))
        fig_diff.savefig("%s_%s_npreddiff.%s"%(args.output, k, args.plot))
        figs.append(fig_map)
        figs.append(fig_diff)


if __name__ == '__main__':
    main()
