#! /usr/bin/env python

"""
Script for generating mock data in Athena++ HDF5 (.athdf) format.
"""

# Python standard modules
import argparse

# Numerical modules
import numpy as np

# Other modules
import h5py

# Main function
def main(**kwargs):

  # Construct grid
  rf = np.exp(np.linspace(np.log(kwargs['r_min']), np.log(kwargs['r_max']), kwargs['n_r'] + 1))
  thf = np.linspace(0.0, np.pi, kwargs['n_th'] + 1)
  phf = np.linspace(0.0, 2.0 * np.pi, kwargs['n_ph'] + 1)
  r = 0.5 * (rf[:-1] + rf[1:])
  th = 0.5 * (thf[:-1] + thf[1:])
  ph = 0.5 * (phf[:-1] + phf[1:])

  # Calculate cutoff
  cutoff_r = np.where((r < kwargs['cutoff_r_min']) | (r > kwargs['cutoff_r_max']), 0.0, 1.0)
  cutoff_th = \
      np.where((th < kwargs['cutoff_th_min']) | (th > np.pi - kwargs['cutoff_th_min']), 0.0, 1.0)
  cutoff_ph = np.ones_like(ph)
  cutoff = cutoff_r[None,None,:] * cutoff_th[None,:,None] * cutoff_ph[:,None,None]

  # Calculate perturbation
  pert_r = np.cos(2.0 * np.pi * kwargs['pert_n_r'] * np.log(r / kwargs['cutoff_r_min']) \
      / np.log(kwargs['cutoff_r_max'] / kwargs['cutoff_r_min']))
  pert_th = -np.cos(2.0 * np.pi * kwargs['pert_n_th'] * (th - kwargs['cutoff_th_min']) \
      / (np.pi - 2.0 * kwargs['cutoff_th_min']))
  pert_ph = np.cos(kwargs['pert_n_ph'] * ph)
  pert = \
      1.0 + kwargs['pert_amp'] * pert_r[None,None,:] * pert_th[None,:,None] * pert_ph[:,None,None]

  # Calculate cell values
  rho = kwargs['rho_amp'] * r[None,None,:] ** -kwargs['rho_r_power'] \
      * np.exp(-np.abs(th[None,:,None] - np.pi / 2.0) / kwargs['rho_th_scale']) * pert * cutoff
  rho = np.maximum(rho, kwargs['rho_floor'])
  pgas = kwargs['pgas_amp'] * r[None,None,:] ** -kwargs['pgas_r_power'] \
      * np.exp(-np.abs(th[None,:,None] - np.pi / 2.0) / kwargs['pgas_th_scale']) * pert ** 2 \
      * cutoff
  pgas = np.maximum(pgas, kwargs['pgas_floor'])
  uur = np.zeros_like(rho)
  uuth = np.zeros_like(rho)
  uuph = kwargs['uph_amp'] * r[None,None,:] ** -kwargs['uph_r_power'] \
      * np.exp(-np.abs(th[None,:,None] - np.pi / 2.0) / kwargs['uph_th_scale']) * cutoff
  rr = np.maximum(r[None,None,:] * np.sin(th[None,:,None]), kwargs['cutoff_r_min'])
  bbz = kwargs['Bz_amp'] * rr ** -kwargs['Bz_R_power']
  bbr = np.cos(th[None,:,None]) * bbz * np.ones_like(ph[:,None,None])
  bbth = -np.sin(th[None,:,None]) / r[None,None,:] * bbz * np.ones_like(ph[:,None,None])
  bbph = kwargs['Bph_amp'] * r[None,None,:] ** -kwargs['Bph_r_power'] \
      * np.exp(-np.abs(th[None,:,None] - np.pi / 2.0) / kwargs['Bph_th_scale']) \
      * np.ones_like(ph[:,None,None])
  if not kwargs['Bph_no_flip']:
    bbph *= np.where(th > np.pi / 2.0, -1.0, 1.0)[None,:,None]

  # Open file for writing
  with h5py.File(kwargs['filename'], 'w') as f:

    # Write file-level attributes
    f.attrs.create('NumCycles', 0, dtype=np.int32)
    f.attrs.create('Time', 0.0, dtype=np.float32)
    f.attrs.create('Coordinates', 'kerr-schild', dtype='|S11')
    f.attrs.create('RootGridX1', (rf[0], rf[-1], (rf[-1] / rf[0]) ** (1.0 / len(r))),
        dtype=np.float32)
    f.attrs.create('RootGridX2', (thf[0], thf[-1], 1.0), dtype=np.float32)
    f.attrs.create('RootGridX3', (phf[0], phf[-1], 1.0), dtype=np.float32)
    f.attrs.create('RootGridSize', (len(r), len(th), len(ph)), dtype=np.int32)
    f.attrs.create('NumMeshBlocks', 1, dtype=np.int32)
    f.attrs.create('MeshBlockSize', (len(r), len(th), len(ph)), dtype=np.int32)
    f.attrs.create('MaxLevel', 0, dtype=np.int32)
    f.attrs.create('NumVariables', [5, 3], dtype=np.int32)
    f.attrs.create('DatasetNames', ['prim', 'B'], dtype='|S21')
    f.attrs.create('VariableNames',
        ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3'], dtype='|S21')

    # Write grid data
    f.create_dataset('Levels', data=(0,), dtype=np.int32)
    f.create_dataset('LogicalLocations', data=(0,0,0), shape=(1,3), dtype=np.int64)
    f.create_dataset('x1f', data=rf, shape=(1,len(rf)), dtype=np.float32)
    f.create_dataset('x2f', data=thf, shape=(1,len(thf)), dtype=np.float32)
    f.create_dataset('x3f', data=phf, shape=(1,len(phf)), dtype=np.float32)
    f.create_dataset('x1v', data=r, shape=(1,len(r)), dtype=np.float32)
    f.create_dataset('x2v', data=th, shape=(1,len(th)), dtype=np.float32)
    f.create_dataset('x3v', data=ph, shape=(1,len(ph)), dtype=np.float32)

    # Write cell data
    f.create_dataset('prim', shape=(5,1,len(ph),len(th),len(r)), dtype=np.float32)
    f['prim'][0,0,...] = rho
    f['prim'][1,0,...] = pgas
    f['prim'][2,0,...] = uur
    f['prim'][3,0,...] = uuth
    f['prim'][4,0,...] = uuph
    f.create_dataset('B', shape=(3,1,len(ph),len(th),len(r)), dtype=np.float32)
    f['B'][0,0,...] = bbr
    f['B'][1,0,...] = bbth
    f['B'][2,0,...] = bbph

# Parse inputs and execute main function
if __name__ == '__main__':

  # Prepare parser
  parser = argparse.ArgumentParser()

  # Prepare for filename input
  parser.add_argument('filename', help='name of simulation data file to write')

  # Prepare for grid inputs
  r_min = 2.0 * 25.0 ** (-1.0 / 75.0)
  r_max = 2.0 * 25.0 ** (76.0 / 75.0)
  parser.add_argument('--r_min', type=float, default=r_min, help='minimum radial coordinate')
  parser.add_argument('--r_max', type=float, default=r_max, help='maximum radial coordinate')
  parser.add_argument('--n_r', type=int, default=77, help='number of cells in radial direction')
  parser.add_argument('--n_th', type=int, default=64, help='number of cells in polar direction')
  parser.add_argument('--n_ph', type=int, default=128,
      help='number of cells in azimuthal direction')

  # Prepare for density inputs
  parser.add_argument('--rho_amp', type=float, default=1.0, help='density coefficient')
  parser.add_argument('--rho_r_power', type=float, default=0.5,
      help='density power-law index in radius')
  parser.add_argument('--rho_th_scale', type=float, default=np.pi/8.0,
      help='density scale height in polar angle')
  parser.add_argument('--rho_floor', type=float, default=1.0e-8, help='density floor')

  # Prepare for gas pressure inputs
  parser.add_argument('--pgas_amp', type=float, default=0.1, help='gas pressure coefficient')
  parser.add_argument('--pgas_r_power', type=float, default=1.5,
      help='gas pressure power-law index in radius')
  parser.add_argument('--pgas_th_scale', type=float, default=np.pi/8.0,
      help='gas pressure scale height in polar angle')
  parser.add_argument('--pgas_floor', type=float, default=1.0e-9, help='gas pressure floor')

  # Prepare for azimuthal velocity inputs
  r_isco = 6.0
  omega_isco = r_isco ** -1.5
  gamma_isco = (1.0 - 2.0 / r_isco - r_isco ** 2 * omega_isco ** 2) ** -0.5
  uph_r_power = 1.5
  uph_amp = gamma_isco * omega_isco * r_isco ** uph_r_power
  parser.add_argument('--uph_amp', type=float, default=uph_amp,
      help='azimuthal velocity coefficient')
  parser.add_argument('--uph_r_power', type=float, default=uph_r_power,
      help='azimuthal velocity power-law index in radius')
  parser.add_argument('--uph_th_scale', type=float, default=np.pi/8.0,
      help='azimuthal velocity scale height in polar angle')

  # Prepare for azimuthal magnetic field inputs
  parser.add_argument('--Bph_amp', type=float, default=0.2,
      help='azimuthal magnetic field coefficient')
  parser.add_argument('--Bph_r_power', type=float, default=1.75,
      help='azimuthal magnetic field power-law index in radius')
  parser.add_argument('--Bph_th_scale', type=float, default=np.pi/8.0,
      help='azimuthal magnetic field scale height in polar angle')
  parser.add_argument('--Bph_no_flip', action='store_true',
      help='do not flip azimuthal magnetic field across midplane')

  # Prepare for vertical magnetic field inputs
  parser.add_argument('--Bz_amp', type=float, default=0.02,
      help='vertical magnetic field coefficient')
  parser.add_argument('--Bz_R_power', type=float, default=0.75,
      help='vertical magnetic field power-law index in cylindrical radius')

  # Prepare for cutoff inputs
  parser.add_argument('--cutoff_r_min', type=float, default=2.0, help='inner cutoff')
  parser.add_argument('--cutoff_r_max', type=float, default=50.0, help='outer cutoff')
  parser.add_argument('--cutoff_th_min', type=float, default=np.pi/16.0,
      help='polar angle cutoff')

  # Prepare for perturbation inputs
  pert_r_num_help = \
      'number (possibly fractional) of perturbation wavelengths in (logarithmic) radial coordinate'
  parser.add_argument('--pert_amp', type=float, default=0.1,
      help='amplitude of density and temperature perturbation')
  parser.add_argument('--pert_n_r', type=float, default=3.0, help=pert_r_num_help)
  parser.add_argument('--pert_n_th', type=float, default=2.0,
      help='number (possibly fractional) of perturbation wavelengths in polar coordinate')
  parser.add_argument('--pert_n_ph', type=int, default=4,
      help='integer number of perturbation wavelengths in azimuthal coordinate')

  # Parse inputs
  args = parser.parse_args()

  # Execute main function
  main(**vars(args))
