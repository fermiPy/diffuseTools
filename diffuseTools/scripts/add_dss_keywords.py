#!/usr/bin/env python
#
import os
import sys

import argparse

from astropy.io import fits

USAGE = "add_dss_keywords.py [options]" 
DESCRIPTION = "Copy DSS keywords from a counts map to an output file"

PARSER = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
PARSER.add_argument('--cmap', type=str,default=None,
                    help='Input counts cube file')
PARSER.add_argument('--outfile', type=str,default=None,
                    help='Output file')

def copy_dss_keywords(hin, hout):
    for k in hin.keys():
        if k == 'NDSKEYS':
            hout[k] = hin[k]
        elif k[0:2] == "DS":
            hout[k] = hin[k]

def main():
    args = PARSER.parse_args(sys.argv[1:])
    fcmap = fits.open(args.cmap)
    fout = fits.open(args.outfile, 'update')
    
    print ("Copy DSS keywords for PRIMARY from %s to %s"%(args.cmap, args.outfile))
    copy_dss_keywords(fcmap[0].header, fout[0].header)
    if fcmap[1].name == 'SKYMAP':
        print ("Copy DSS keywords for SKYMAP from %s to %s"%(args.cmap, args.outfile))
        copy_dss_keywords(fcmap[1].header, fout[1].header)
    fout.flush()

if __name__ == '__main__':
    main()

