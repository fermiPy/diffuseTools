#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

import sys

import argparse
from fermipy import utils
from diffuseTools.plot_tools import plot_map_sums, plot_map_diff

def main():
    usage = "make_spectra.py [options]" 
    description = "Plot spectral components from a results file"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('-i', '--input', type=str, default=None,
                        help='Input file')    
    parser.add_argument('-c', '--compare', type=str, default=None,
                        help='Input file to compare')    
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output file prefix')
    parser.add_argument('-p', '--plot', type=str, default='png',
                        help='Plot type')
    parser.add_argument('-d', '--comp_dict', type=str, default=None,
                        help='Component dict')

    args = parser.parse_args(sys.argv[1:])

    if args.input is not None:
        map_dict = utils.load_yaml(args.input)
    else:
        map_dict = None

    if args.compare is not None:
        compare_map_dict = utils.load_yaml(args.compare)
    else:
        compare_map_dict = None
    
    comp_dict = utils.load_yaml(args.comp_dict)
        
    figs = []

    for k, v in comp_dict.items():
        fig_map = plot_map_sums(map_dict, v, compare_map_dict)
        fig_diff = plot_map_diff(map_dict, v, compare_map_dict)
        fig_map.savefig("%s_%s_mapsum.%s"%(args.output, k, args.plot))
        fig_diff.savefig("%s_%s_mapsumdiff.%s"%(args.output, k, args.plot))
        figs.append(fig_map)
        figs.append(fig_diff)
    
    return figs

if __name__ == '__main__':
    figs = main()
