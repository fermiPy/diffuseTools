#!/usr/bin/env python
#

import argparse
import time
import sys

import numpy as np
from fermipy.gtanalysis import GTAnalysis

from diffuseTools.fermipy_tools import build_srcdict, get_src_names, set_wts_get_npred_wt


def main():

    # Argument defintion
    usage = "usage: %(prog)s [options]" 
    description = "Load a fermipy configuration and snapshot"

    parser = argparse.ArgumentParser(usage,description=description)
    parser.add_argument('-c', "--config",
                        type=argparse.FileType('r'),
                        default="config.yaml",
                        help="Config file")
    parser.add_argument('-i', "--input",
                        type=argparse.FileType('r'),
                        default="baseline.npy",
                        help="Input file")
    parser.add_argument('-o', "--output",
                        type=str,
                        default=None,
                        help="Output file")
    parser.add_argument('-p', "--pars",
                        type=str,
                        default=None,
                        help="Parameter file")
    parser.add_argument('-m', "--mask",
                        type=str,
                        default=None,
                        help="Mask file")

    args = parser.parse_args(sys.argv[1:])

    gta = GTAnalysis.create(args.input.name, args.config.name, 
                            params=args.pars, mask=args.mask)

    for name in gta.like.sourceNames():
        print('Initializing source %s' % name)
        gta._init_source(name)
    gta._update_roi()

    gta.write_roi(args.output, save_model_map=True, save_weight_map=True)

    return gta


if __name__ == '__main__':
    gta = main()
