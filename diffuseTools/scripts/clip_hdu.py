#!/usr/bin/env python

import sys

import argparse
from astropy.io import fits
import numpy as np

def clip_hdu(hdu, minval, maxval):
    """
    """
    for icol, col in enumerate(hdu.columns):
        load = col.array[0]
        hdu.columns[icol].array = col.array.clip(min=minval, max=maxval)
    hdu_out = fits.BinTableHDU.from_columns(hdu.columns, header=hdu.header, name=hdu.name)
    return hdu_out


if __name__ == '__main__':

    usage = "clib_hdu.py [options]" 
    description = "CLIP values in an HDU in a set of files"

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('--hdu', type=str, default=None,
                        help='HDU name')
    parser.add_argument('--min', type=float, default=None,
                        help='Minimum value')
    parser.add_argument('--max', type=float, default=None,
                        help='Maximum value')
    parser.add_argument('files', nargs='+', help='Names of input files')

    args = parser.parse_args(sys.argv[1:])

    for fin in args.files:
        print ("Working on %s"%fin)
        hdulist = fits.open(fin,'update')
        hdu = hdulist[args.hdu]
        hdu_out = clip_hdu(hdu, args.min, args.max)
        hdulist[args.hdu] = hdu_out
        hdulist.flush()
        hdulist.close()

        
