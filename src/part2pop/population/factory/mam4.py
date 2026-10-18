#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a population from a PARTMC NetCDF file
@author: Laura Fierce
"""

from ..base import ParticlePopulation
from part2pop import make_particle_from_masses
from part2pop.species.registry import get_species

import numpy as np
import os
from pathlib import Path
from .registry import register
from ...constants import MOLAR_MASS_DRY_AIR, R, DENSITY_LIQUID_WATER
from ...utilities import power_moments_from_lognormal
from .binned_lognormals import build as build_binned_lognormals


def _select_active_mode_values(values, indices, n_modes, name, allow_scalar=False):
    """Select by original mode index; scalar bin counts apply to every mode."""
    array = np.asarray(values)
    if allow_scalar and array.ndim == 0:
        return array.item()
    if array.ndim != 1 or len(array) != n_modes:
        raise ValueError(f"{name} must have one entry per original MAM mode ({n_modes}).")
    return array[indices].tolist()


try:
    import netCDF4
    _HAS_NETCDF4 = True
except ModuleNotFoundError:
    netCDF4 = None
    _HAS_NETCDF4 = False
if _HAS_NETCDF4:
    @register("mam4") # only registers if netCDF4 is available
    def build(config):
        mam4_dir = Path(config['mam4_dir'])
        timestep = config['timestep']
        GSDs = config['GSD']
        if config.get("D_is_wet", True):
            mode_species = ['SO4', 'OC', 'H2O']
        else:
            mode_species = ['SO4', 'OC']
        aero_spec_names = [mode_species.copy() for _ in GSDs]
        
        N_bins = config['N_bins']
        

        if timestep == 0:
            raise ValueError('timestep=0 is invalid. Specify timestep = 1 for initial conditions')
        elif timestep == 1:
            mam_input = mam4_dir / 'namelist'
            Ns = np.zeros([len(GSDs)])
            
            for kk in range(len(Ns)):
                Ns[kk] = get_mam_input(
                    'numc' + str(kk+1),
                    mam_input)
            
            # GMDs = np.zeros([len(GSDs)])
            # try:
            #     for kk in range(len(GMDs)):
            #         GMDs[kk] = get_mam_input(
            #             'dgn' + str(kk+1),
            #             mam_input)
            # except ValueError:
                # if dgn not in namelist, use default initial GMDs
            GMDs = config.get('GMD_init', [1.1e-7, 2.6e-8, 2e-6, 5e-8]) # default initial GMDs for MAM4 (check fresh BC mode)
            
            aero_spec_fracs = []
            for kk in range(len(Ns)):
                spec_frac = []
                total_frac = 0.
                for spec_name in aero_spec_names[kk]:
                    # fixme: come up with more elegant way to handle PartMC --> MAM4 species name mapping
                    if spec_name.lower() == 'oc':
                        spec_name = 'pom'
                    
                    if spec_name.lower() == 'h2o':
                        frac = 0. 
                        #raise Warning('H2O mass fraction cannot be specified in MAM4 namelist; it is computed internally based on other species and hygroscopic growth. Setting frac_H2O=0 here.')
                    else: 
                        frac = get_mam_input(
                            'mf' + spec_name.lower() + str(kk+1),
                            mam_input)
                                        
                    spec_frac.append(frac)
                    total_frac += frac
                spec_frac = np.array(spec_frac)
                spec_frac /= total_frac
                aero_spec_fracs.append(spec_frac)
            aero_spec_fracs = np.array(aero_spec_fracs)

            lognormals_cfg = {
                'type': 'binned_lognormals',
                'N_sigmas': config.get('N_sigmas',6),
                'D_min':config.get('D_min',None), #fixme: option for +/- sigmas
                'D_max':config.get('D_max',None),
                'N_bins': N_bins,
                'N': Ns,
                'GMD': GMDs,
                'GSD': GSDs,
                'aero_spec_names': aero_spec_names,
                'aero_spec_fracs': aero_spec_fracs,
                }
        else:
            output_tt_idx = timestep - 2 # timestep=2 is output index 0 = after 1 time step

            # fixme: make this +/- certain number of sigmas (rather than min/max diams)
            
            N_bins = config['N_bins']
            output_filename = mam4_dir / 'mam_output.nc'
            with netCDF4.Dataset(output_filename) as currnc:
                num_aer = currnc['num_aer'][:, output_tt_idx]
                so4_aer = currnc['so4_aer'][:, output_tt_idx]
                soa_aer = currnc['soa_aer'][:, output_tt_idx]
                dgn_a = currnc['dgn_a'][:, output_tt_idx]
                dgn_awet = currnc['dgn_awet'][:, output_tt_idx]

            p = config['p']
            T = config['T']
            #from .binned_lognormals import build as build_binned_lognormals
            
            lognormals_cfg = config

            # fixme: make this right
            rho_dry_air = MOLAR_MASS_DRY_AIR * p / (R * T)
            Ns = num_aer * rho_dry_air # N/m^3
            idx, = np.where(Ns>0)
            # Every per-mode input refers to the original file mode order.
            # Apply the same indices before moments or population construction;
            # truncating a prefix would misidentify modes after an inactive one.
            GSDs = _select_active_mode_values(GSDs, idx, len(Ns), 'GSD')
            N_bins = _select_active_mode_values(
                N_bins, idx, len(Ns), 'N_bins', allow_scalar=True)
            Ns = Ns[idx]
            D_is_wet = config.get("D_is_wet", True)
            GMDs_wet = dgn_awet[idx]
            GMDs = dgn_a[idx]
            
            mass_so4 = so4_aer[idx] * rho_dry_air
            mass_soa = soa_aer[idx] * rho_dry_air
            # fixme: not actually sure what assumptions go into this for MAM4, assuming same GSD to compute mass H2O
            mass_h2o = DENSITY_LIQUID_WATER * np.pi/6 * np.array(
                    [
                        (power_moments_from_lognormal(3,N,gmd_wet,gsd) 
                        - power_moments_from_lognormal(3, N, gmd, gsd)) 
                        for (N, gmd, gmd_wet, gsd) in zip(Ns, GMDs, GMDs_wet, GSDs)])
            if D_is_wet:
                aero_spec_names = [mode_species.copy() for _ in idx]
                aero_spec_fracs = []
                for (m_so4,m_soa,m_h2o) in zip(mass_so4,mass_soa,mass_h2o):
                    spec_frac = np.array([m_so4,m_soa,m_h2o])
                    spec_frac /= np.sum(spec_frac)
                    spec_frac[np.isnan(spec_frac)] = 0.
                    aero_spec_fracs.append(spec_frac)
            else:
                aero_spec_names = [mode_species.copy() for _ in idx]
                aero_spec_fracs = []
                for (m_so4,m_soa) in zip(mass_so4,mass_soa):
                    spec_frac = np.array([m_so4,m_soa])
                    spec_frac /= np.sum(spec_frac)
                    spec_frac[np.isnan(spec_frac)] = 0.
                    aero_spec_fracs.append(spec_frac)

            if D_is_wet:
                GMDs_cfg = GMDs_wet
            else:
                GMDs_cfg = GMDs
            lognormals_cfg = {
                "type": "binned_lognormals",
                "N": Ns,
                "N_sigmas": config.get('N_sigmas', 6),
                "D_is_wet": D_is_wet,
                "GMD": GMDs_cfg,
                "GSD": GSDs,
                "aero_spec_names": aero_spec_names,
                "aero_spec_fracs": aero_spec_fracs,
                "N_bins": N_bins,
                "D_min": config.get('D_min', None),
                "D_max": config.get('D_max', None),
            }
        # Preserve explicit species properties and the overrides inherited by
        # downstream optics/freezing, as the other population factories do.
        lognormals_cfg['species_modifications'] = config.get('species_modifications', {})
        mam4_population = build_binned_lognormals(lognormals_cfg)
        return mam4_population

else:
    def build(config):
        raise ModuleNotFoundError(
            "Install netCDF4 to read MAM4 files: generate the environment-partmc.yml file using tools/create_conda_env.py, "
            "then create and activate the 'part2pop' conda environment (conda env create -f environment-partmc.yml -n part2pop)."
        )

# fimxe: better way to deal with MAM input? send in namelist (or path to namelist)
# import f90nml
# namelist = f90nml.read('namelist')
# qso2 = namelist['chem_input']['qso2']
def get_mam_input(varname,mam_input_filename):

    f_input = open(mam_input_filename,'r')
    input_lines = f_input.readlines()
    yep = 0
    for oneline in input_lines:
        if varname in oneline:
            yep += 1 
            idx1,=np.where([hi == '=' for hi in oneline])
            idx2,=np.where([hi == ',' for hi in oneline])
            vardat = float(''.join([oneline[ii] for ii in range(idx1[0]+1,idx2[0])]))
    if yep == 0:
        # raise ValueError(varname,'is not a MAM input parameter')
        vardat = 0.0
        print(varname,'is not a MAM input parameter; returning 0.0')    
    elif yep > 1:
        raise ValueError('more than one line in ', mam_input_filename, 'starts with', varname)
    return vardat