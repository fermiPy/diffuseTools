#!/usr/bin/env python
#

import argparse
import time
import sys

import numpy as np
from fermipy.gtanalysis import GTAnalysis

def main():
    # Argument defintion
    usage = "usage: %(prog)s [options]" 
    description = "Run an analysis script"
    
    parser = argparse.ArgumentParser(usage,description=description)
    parser.add_argument('-c', "--config",
                        type=argparse.FileType('r'),
                        default="config.yaml",
                        help="Input file")
    parser.add_argument('-o', "--output",
                        type=str,
                        default='baseline',
                        help="Output file prefix")

    args = parser.parse_args(sys.argv[1:])
   
    gta = GTAnalysis(args.config.name)
    gta.setup()
    
    gta.write_roi(args.output, save_model_map=True, save_weight_map=True, make_plots=True)
    return gta


if __name__ == '__main__':
    gta = main()
