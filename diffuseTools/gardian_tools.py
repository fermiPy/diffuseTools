# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tools to interact with GARDIAN
"""
from __future__ import absolute_import, division, print_function

import sys

import argparse
import yaml
import numpy as np
from scipy.interpolate import interp1d
from fermipy.skymap import read_map_from_fits


def spec_type(srcname, conv_dict):
    if srcname[0:4] == 'FL8Y':
        return 'FL8Y'
    elif conv_dict.has_key(srcname):
        return conv_dict[srcname]['spec_type']
    raise ValueError("Can't match %s" % srcname)
                                      

def read_lines(filepath):    
    with open(filepath, 'r') as f:
        l = [line.strip() for line in f]
    return l

def get_status(lines):
    if lines[0].find('successful') > 0:
        success = True
    else:
        success = False
    mll = float(lines[-2].split()[4])
    return success, mll

def read_token(token):
    if token == 'N/A':
        return np.nan
    return float(token)

def read_prior(token):
    if token == 'N/A':
        return None
    return token

def read_val_and_error(token):
    sub_tokens = token.split('(')
    value = float(sub_tokens[0].strip())
    if len(sub_tokens) == 1:
        error = np.nan
    else:
        error = float(sub_tokens[1][0:-1].strip())
    return value, error

def read_var_line(l):
    tokens = l.split('|')
    if len(tokens) != 7:
        print("Wrong number of tokens, %i %s"%(len(tokens), l))
    var_name = tokens[0].strip()
    value, error = read_val_and_error(tokens[1].strip())
    min_val = read_token(tokens[2].strip())
    max_val = read_token(tokens[3].strip())
    prior = read_prior(tokens[4].strip())
    par1 = read_token(tokens[5].strip())
    par2 = read_token(tokens[6].strip())
    d = dict(value=value, error=error,
             min=min_val, max=max_val, 
             prior=prior,
             par1=par1, par2=par2)
    return var_name, d

def read_comp_names(line):
    return line.split()[3:]

def read_values(line):
    return np.array(read_comp_names(line)).astype(float)

def read_energies(line):
    return np.array(line.split()[0:2]).astype(float)

def read_gardian_counts(filepath):
    lines = read_lines(filepath)
    comp_names = read_comp_names(lines[1])
    vals = []
    evals = []
    for l in lines[2:]:
        vals += [read_values(l)]
        evals += [read_energies(l)]
    vv = np.vstack(vals).T
    ee = np.vstack(evals).T
    energies = np.sqrt(ee[0]*ee[1])
    o = dict(Energies=energies)
    for comp_name, val in zip(comp_names, vv):
        o[comp_name] = val
    return o    

def make_npred_dict(inputPrefix, conv_dict):
    o = {}
    comps = conv_dict['COMPONENTS']
    src_dict = conv_dict['SRC_NAME_DICT']
    for k in comps:
        fname_gard = "%s_ModelComponentsCounts/%s.txt"%(inputPrefix, k)
        counts = read_gardian_counts(fname_gard)
        for kk,vv in sorted(counts.items()):
            if kk=="Energies":
                src = kk
            else:
                src = src_dict[kk]
            if not o.has_key(src):
                o[src] = []
            o[src] += [vv]
    return o

def get_var_dict(lines):
    d = {}
    for l in lines[4:-3]:
        name, dd = read_var_line(l)
        d[name] = dd    
    return d

def get_src_and_varnames(k):
    i = k.rfind('_')
    src_name = k[0:i]
    var_name = k[i+1:]
    if k[0:3] == 'iso':
        src_name = 'iso_map'
        var_name = "%s_%s"%(k.split('_')[1],var_name)
    elif k[0:10] == 'Scalepatch':
            src_name = 'patches'
            var_name = "%s_%s"%(var_name, k[10])
    return (src_name, var_name)

def get_src_dict(var_dict):
    d = {}
    for k,v in var_dict.items():
        src_name, var_name = get_src_and_varnames(k)
        if src_name == 'fixed_negative':
            continue
        if not d.has_key(src_name):
            d[src_name] = {}
        d[src_name][var_name] = v
    return d

def translate_lp(var_dict):
    Eb_dict = dict(value=5.62341325e+03)
    beta_def = dict(value=0.0)
    out_dict = dict(norm=var_dict['pref'],
                    alpha=var_dict['index'],
                    beta=var_dict.get('beta', beta_def),
                    Eb=Eb_dict)
    out_dict['beta']['value'] *= -1.
    out_dict['alpha']['value'] *= -1.    
    return out_dict

def translate_pl(var_dict):
    Scale_dict = dict(value=5.477225e+03)
    out_dict = dict(Prefactor=var_dict['pref'],
                    Index=var_dict['index'],
                    Scale=Scale_dict)
    out_dict['Index']['value'] *= -1.
    return out_dict

def translate_constant(var_dict):
    out_dict = dict(Value=var_dict['pref'])
    return out_dict

def translate_fixed_srcs(var_dict):
    out_dict = dict(Value=var_dict['pref'])
    return out_dict

def scale_bands(var_dict, n_band):
    out_dict = {}
    norm = var_dict['dNdE0']['value']    
    for i in range(n_band):
        out_dict['dNdE%i' % i] = var_dict['dNdE%i' % i]
        if i != 0:
            out_dict['dNdE%i' % i]['value'] /= norm
    return out_dict

def translate_bands(var_dict, n_band):
    out_dict = {}
    for i in range(n_band):
        try:
            out_dict['dNdE%i' % i] = var_dict['int%02i' % i]
        except KeyError:
            print (var_dict)
            return var_dict
    out_dict = scale_bands(out_dict, n_band)
    return out_dict


def translate_patches(var_dict):
    out_dict = dict(dNdE0=var_dict['pref_0'],
                    dNdE1=var_dict['pref_1'],
                    dNdE2=var_dict['pref_2'],
                    dNdE3=var_dict['pref_3'])
    out_dict = scale_bands(out_dict, 4)
    return out_dict

def translate_iso(var_dict):
    out_dict = dict(dNdE0=var_dict['PSF3_int00'],
                    dNdE1=var_dict['PSF3_int01'],
                    dNdE2=var_dict['PSF23_int00'],
                    dNdE3=var_dict['PSF23_int01'],
                    dNdE4=var_dict['PSF123_int00'],
                    dNdE5=var_dict['PSF123_int01'],
                    dNdE6=var_dict['PSF0123_int00'],
                    dNdE7=var_dict['PSF0123_int01'],
                    dNdE8=var_dict['PSF0123_int02'],
                    dNdE9=var_dict['PSF0123_int03'],
                    dNdE10=var_dict['PSF0123_int04'],
                    dNdE11=var_dict['PSF0123_int05'],
                    dNdE12=var_dict['PSF0123_int06'],
                    dNdE13=var_dict['PSF0123_int07'],
                    dNdE14=var_dict['PSF0123_int08'])
    out_dict = scale_bands(out_dict, 15)
    return out_dict

def translate_FL8Y(var_dict):
    out_dict = dict(pref=var_dict['pref'])
    return out_dict

def translate_var_dict(var_dict, spectype):
    if spectype=="LogParabola":
        return translate_lp(var_dict)
    elif spectype=="ConstantValue":
        return translate_constant(var_dict)
    elif spectype=="Powerlaw":
        return translate_pl(var_dict)
    elif spectype=="BinByBin5":
        return translate_bands(var_dict, 5)
    elif spectype=="BinByBin9":
        return translate_bands(var_dict, 9)
    elif spectype=="Iso":
        return translate_iso(var_dict)
    elif spectype=="Patches":
        return translate_patches(var_dict)
    elif spectype=="FL8Y":
        return translate_FL8Y(var_dict)
    raise ValueError("Did not recognize spectrum %s" % spectype)


def translate_src_dict(src_dict, conv_dict):
    out_dict = {}
    for k, v in sorted(src_dict.items()):
        if k[0:4] == 'FL8Y': 
            src = k
        elif k == 'fixed_negative':
            src = k
        else:
            src = conv_dict[k]['name']
        if src in [None, 'none', 'None']:
            continue
        spectype = spec_type(k, conv_dict)
        try:
            out_dict[src] = translate_var_dict(v, spectype)
        except ValueError as msg:
            print (spectype, k, v)
            raise ValueError(msg)
        try:
            out_dict[src]['SpectrumType'] = spectype
        except TypeError:
            print (out_dict[src], k, v)
    return out_dict

def adjust_xco_norm(src_dict, k, kk, v, vv):
    CO_pref_dict = v.pop('Value')
    HI_pref_dict = vv['Prefactor']
    HI_val = HI_pref_dict['value']
    CO_val = CO_pref_dict['value']
    CO_err = CO_pref_dict['error']
    src_dict[k]['SpectrumType'] = src_dict[kk]['SpectrumType'] 
    src_dict[k]['Index'] = src_dict[kk]['Index'].copy()
    src_dict[k]['Prefactor'] = src_dict[kk]['Prefactor'].copy()
    src_dict[k]['Scale'] = src_dict[kk]['Scale'].copy()
    print ("Adjusting XCO", k, kk, HI_val, CO_val)
    src_dict[k]['Prefactor']['value'] *= CO_val

def adjust_xco_bands(src_dict, k, kk, v, vv):
    CO_pref_dict = v.pop('Value')
    CO_val = CO_pref_dict['value']
    CO_err = CO_pref_dict['error']
    for kc, vc in src_dict[kk].items():
        if isinstance(vc, str):
            src_dict[k][kc] = vc
        elif isinstance(vc, dict):
            src_dict[k][kc] = vc.copy()
    var = 'dNdE0'
    var_dict = vv[var]
    HI_val = var_dict['value']
    print ("Adjusting XCO", k, kk, var, HI_val, CO_val)
    src_dict[k][var]['value'] *= CO_val
        
    #print ("Adjusting XCO", k, kk, HI_val, CO_val)
    #src_dict[k]['Prefactor']['value'] *= CO_val


def adjust_xco(src_dict):    
    for k, v in sorted(src_dict.items()):
        if k.find('merged_CO') < 0:
            continue
        kk = k.replace('CO','HI')
        vv = src_dict[kk]
        if vv.has_key('Prefactor'):
            adjust_xco_norm(src_dict, k, kk, v, vv)
        else:
            adjust_xco_bands(src_dict, k, kk, v, vv)


def adjust_fl8y(src_dict, cat_dict):
    if cat_dict is None:
        return src_dict
    out_dict = cat_dict.copy()
    for k,v in src_dict.items():
        if out_dict.has_key(k):
            st = out_dict[k]['SpectrumType']
            if st == 'LogParabola':
                norm_var = 'norm'
            else:
                norm_var = 'Prefactor'
            try:
                print ("Scale ", k, out_dict[k][norm_var]['value'],  v['pref']['value'])
                out_dict[k][norm_var]['value'] *= v['pref']['value']
            except KeyError as msg:
                print (k,v,out_dict[k])
                raise KeyError(msg)
        else:
            out_dict[k] = v
    return out_dict

def opt_write_val(key, d, stream, prefix="    "):
    if d.has_key(key) and np.isfinite(d[key]):
        stream.write("%s%s : %.4e\n"%(prefix, key, d[key]))

def opt_write_val_short(key, d, stream):
    if d.has_key(key) and np.isfinite(d[key]):
        stream.write("%0.4e"%(d[key]))
    else:
        stream.write('nan')


def pretty_print_input(src_dict, outpath):
    if outpath is None:
        out = sys.stdout
    else:
        out = open(outpath, 'w')

    for src, var_dict in sorted(src_dict.items()):
        out.write("%s: \n"%(src))
        spectype = var_dict.pop('SpectrumType')
        out.write('  SpectrumType : %s\n'%(spectype))
        for var_name, var_vals in sorted(var_dict.items()):
            out.write("  %s:\n"%var_name)
            opt_write_val('value', var_vals, out)
            opt_write_val('error', var_vals, out)
            opt_write_val('min', var_vals, out)
            opt_write_val('max', var_vals, out)
            prior = var_vals.get('prior')
            par1 = var_vals.get('par1')
            par2 = var_vals.get('par2')
            if prior is None:
                continue
            #out.write("    prior : %s(%s, %s)\n"%(prior, par1, par2))
                 

def pretty_print_result(src_dict, outpath):
    if outpath is None:
        out = sys.stdout
    else:
        out = open(outpath, 'w')

    for src, var_dict in sorted(src_dict.items()):
        out.write("%s:\n"%(src))
        spectype = var_dict.pop('SpectrumType')
        out.write('  SpectrumType : %s\n'%(spectype))
        for var_name, var_vals in sorted(var_dict.items()):
            out.write("  %s: ["%(var_name))
            opt_write_val_short('value', var_vals, out)
            out.write(', ')
            opt_write_val_short('error', var_vals, out)
            out.write(', 1.0000e+00, ')
            opt_write_val_short('min', var_vals, out)
            out.write(', ')            
            opt_write_val_short('max', var_vals, out)
            out.write(']\n')
             

def interpolate_factors(bin_energies,ppl_factors):
    PPL_EDGES = np.logspace(np.log10(30), 6, 6)
    PPL_CENTERS = np.sqrt(PPL_EDGES[0:-1]*PPL_EDGES[1:])
    interp = interp1d(PPL_CENTERS, ppl_factors, fill_value='extrapolate')
    l_out = [ interp(bin_evals) for bin_evals in bin_energies ]
    return l_out

def combine_dnm_maps(out_dict, use_factors=False):
    try:
        dnm_pos = out_dict.pop('dnm_pos')
    except KeyError as msg:
        print (out_dict.keys())
        raise KeyError(msg)

    try:
        dnm_neg = out_dict.pop('dnm_neg')
    except KeyError as msg:
        print (out_dict.keys())
        raise KeyError(msg)
       
    bin_energes = out_dict['Energies']
    dnm_merged = []

    if use_factors:
        fp = np.array([0.27647, 0.48268, 0.69477, 0.42879, 0.18468])
        fn = np.array([0.39485, 0.52235, 0.80803, 0.98369, 1.2549])
        fp_a = interpolate_factors(bin_energes, fp)
        fn_a = interpolate_factors(bin_energes, fn)
    else:
        fp_a = np.ones((len(bin_energes)))
        fn_a = np.ones((len(bin_energes)))

    for dp, dn, fpv, fnv in zip(dnm_pos, dnm_neg, fp_a, fn_a):
        dm = (fpv*np.array(dp)) - (fnv*np.array(dn))
        dnm_merged.append(dm)
    out_dict['dnm_merged'] = dnm_merged


def get_map_sums(prefix, bin_comps, model_comp, sel_dep_comps):
    o = []
    for k in bin_comps:
        if model_comp in sel_dep_comps:
            fname_gard = "%s_ModelComponentsCounts/%s_%s%s.fits.gz"%(prefix, k, model_comp, k[0:2])
        else:
            fname_gard = "%s_ModelComponentsCounts/%s_%s.fits.gz"%(prefix, k, model_comp)
        try:
            m = read_map_from_fits(fname_gard)
        except IOError:
            print ("Failed to read %s" % fname_gard)
            o.append([0.])
        o.append(m.counts.sum(1))
    return o

def make_sum_dict(prefix, bin_comps, src_dict, sel_dep_comps):
    o = {}
    sys.stdout.write("Reading %i components\n"%(len(src_dict)))
    for k,v in sorted(src_dict.items()):
        sys.stdout.write('.')
        sys.stdout.flush()
        o[v] = get_map_sums(prefix, bin_comps, k, sel_dep_comps)
    sys.stdout.write('!\n')
    return o

def get_expmap_sums(filepath, energies):
    o = []
    m = read_map_from_fits(filepath)
    exps = m.counts.sum(1)
    vv = np.interp(energies, m.hpx.evals, exps)
    return vv

def read_egb_file(filepath):
    with open(filepath, 'r') as f:
        tokens = [line.split() for line in f]
    v = np.array(tokens).astype(float)
    return v.T[0:2]

def get_int_egb(filepath, energies):
    egb_data = read_egb_file(filepath)
    vv = np.interp(energies, egb_data[0], egb_data[1])
    return vv
