#!/usr/bin/env python
#

# Adapted from 'SummaryPdfFile7.py by Seth Digel and Suttiwat "Bus" Madlee/

import sys
import argparse

import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import healpy as hp

from diffuseTools.plot_tools import make_grid, make_regions, figure_model, figure_frac_resid, figure_signif,\
    figure_region_spectra, figure_region_resid

NSIDE = 32
REGIONS = make_regions()
NCOMP = 4

def main():
    usage = "usage: %(prog)s [options]" 
    description = "Make output summary plots of all-sky fitting"
    
    parser = argparse.ArgumentParser(usage,description=description)
    parser.add_argument('-i', "--input",
                        type=str,
                        default="mcube_baseline",
                        help="Input file prefix")
    parser.add_argument('-c', "--ccube",
                        type=str,
                        default="ccube.fits",
                        help="Input file prefix")
    args = parser.parse_args(sys.argv[1:])

    fontsize = 9
    matplotlib.rcParams.update({'font.size': fontsize})
    
    #pdf_pages = PdfPages("%s.pdf"%(args.input))
    
    ccube = hp.mrdfits(args.ccube, hdu=1)[1:]
    ebins = hp.mrdfits(args.ccube, hdu=2)
    eranges = np.vstack([ebins[1],ebins[2]]).T / 1e6  # keV to GeV
    ecents = np.sqrt(ebins[1]*ebins[2])
    ewidths = ebins[2] - ebins[1]
    
    nebins = len(eranges)
    counts_spectrum = np.zeros((nebins))

    idx = 0
    sys.stdout.write("Making plots: ")

    region_spectra = {}
    
    for i in range(NCOMP):
        mcube_file = "%s_%02i.fits"%(args.input, i)
        mcube = hp.mrdfits(mcube_file, hdu=1)[1:]
       
        for j in range(len(mcube)):
            sys.stdout.write('.')
            sys.stdout.flush()
            
            counts_cast = hp.ud_grade(ccube[idx], NSIDE, power=-2)
            model_cast = hp.ud_grade(mcube[j], NSIDE, power=-2)

            fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
            figure_model(model_cast, eranges[idx])
            figure_frac_resid(counts_cast, model_cast, eranges[idx])

            figure_signif(counts_cast, model_cast, eranges[idx])            
            #pdf_pages.savefig(fig,bbox_inches='tight')
            
            fig.savefig('fig_%s_maps_%02i_%02i.png'%(args.input, i,j), bbox_inches='tight')

            for regName, regPixels in zip(REGIONS['regionNames'], REGIONS['regionPixels']):
                if not region_spectra.has_key(regName):
                    region_spectra[regName] = dict(counts=np.zeros((nebins)),
                                                   model=np.zeros((nebins)))
                region_counts_spectra = region_spectra[regName]['counts']
                region_model_spectra = region_spectra[regName]['model']

                if regPixels is None:
                    region_counts_spectra[idx] = counts_cast.sum()
                    region_model_spectra[idx] = model_cast.sum()
                else:
                    region_counts_spectra[idx] = counts_cast[regPixels].sum()
                    region_model_spectra[idx] = model_cast[regPixels].sum()        

            idx += 1

    fig_region_spectra = figure_region_spectra(REGIONS['regionNames'], region_spectra, ecents, ewidths)
    fig_region_resid = figure_region_resid(REGIONS['regionNames'], region_spectra, ecents, ewidths)

    fig_region_spectra.savefig('fig_%s_region_spectra.png'%(args.input), bbox_inches='tight')
    fig_region_resid.savefig('fig_%s_region_resid.png'%(args.input), bbox_inches='tight')


if __name__ == '__main__':
    main()
