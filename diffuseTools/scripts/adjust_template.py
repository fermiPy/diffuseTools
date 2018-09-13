#!/usr/bin/env python
#

""" 
Adjust a template by multiplying it by a set of scale parametrs
"""

import numpy as np
from fermipy import utils
from gammapy.maps import Map
from scipy.interpolate import interp1d


def make_interpolator(yamlfile):
    scale = utils.load_yaml(yamlfile)
    emids = (np.array(scale['E_min']) + np.array(scale['E_max']))/2.
    func = interp1d(emids, scale['Value'], kind='nearest', bounds_error=False, fill_value='extrapolate')
    return func

def main():

    import os
    import sys
    import argparse

    # Argument defintion
    usage = "usage: %(prog)s [options]" 
    description = "addjust a tempale"

    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument('-o', '--output', type=str, help='Output file', default=None)       
    parser.add_argument('-i', '--input', type=str, help='Input file')
    parser.add_argument('-s', '--scale', type=str, help='File with scaling parameters')
    args = parser.parse_args(sys.argv[1:])

    func = make_interpolator(args.scale)
    input_map = Map.read(args.input)
    energies = input_map.geom.axes[0].center
    
    scale_factors = func(energies)
    output_map = Map.from_geom(input_map.geom)
    output_map.data.flat = (input_map.data.T * scale_factors).T
    
    if args.output is not None:
        print (input_map.geom.conv)
        output_map.write(args.output, conv=input_map.geom.conv)
    else:
        return output_map


if __name__ == "__main__":
    main()
    
