#! /usr/bin/env python

"""
Script for generating mock data to be read by Blacklight.

Supported formats:
  - athena: Standard HDF5 (.athdf) format used by Athena++.
  - iharm3d: HDF5 format with most of the common iharm fields supplied.
  - harm3d: version of older Harm ascii/binary format without excessive duplicate information.
  - harm3d_ext: original version of harm3d, not supported by Blacklight but used by other
      ray-tracing codes.
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
  lr_min = np.log(kwargs['r_min'])
  lr_max = np.log(kwargs['r_max'])
  lrf, dlr = np.linspace(lr_min, lr_max, kwargs['n_r'] + 1, retstep=True)
  rf = np.exp(lrf)
  thf, dth = np.linspace(0.0, np.pi, kwargs['n_th'] + 1, retstep=True)
  phf, dph = np.linspace(0.0, 2.0 * np.pi, kwargs['n_ph'] + 1, retstep=True)
  r = 0.5 * (rf[:-1] + rf[1:])
  th = 0.5 * (thf[:-1] + thf[1:])
  ph = 0.5 * (phf[:-1] + phf[1:])
  lr = np.log(r)
  x2f = thf / np.pi
  x2 = th / np.pi
  dx2 = dth / np.pi

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

  # Calculate extra data for Harm formats
  if kwargs['format'] in ('iharm3d', 'harm3d', 'harm3d_ext'):

    # Calculate metric
    r_hor = 2.0
    sigma = r[None,None,:] ** 2
    f = 2.0 * r[None,None,:] / sigma
    g_tt = -(1.0 - f)
    g_tr = f
    g_tth = 0.0
    g_tph = 0.0
    g_rr = 1.0 + f
    g_rth = 0.0
    g_rph = 0.0
    g_thth = sigma
    g_thph = 0.0
    g_phph = r[None,None,:] ** 2 * np.sin(th[None,:,None]) ** 2
    gtt = -(1.0 + f)
    gtr = f
    gtth = 0.0
    gtph = 0.0
    alpha = 1.0 / np.sqrt(-gtt)
    g00 = gtt
    g01 = gtr / r[None,None,:]
    g02 = gtth / np.pi
    g03 = gtph
    alpha_alt = 1.0 / np.sqrt(-g00)
    det_alt = sigma * np.sin(th[None,:,None]) * r[None,None,:] * np.pi

    # Calculate thermodynamic quantities
    ugas = pgas / (kwargs['gamma_adi'] - 1.0)
    kappa = pgas / rho ** kwargs['gamma_adi']

    # Calculate velocity
    uut = np.sqrt(1.0 + g_rr * uur ** 2 + 2.0 * g_rth * uur * uuth + 2.0 * g_rph * uur * uuph \
        + g_thth * uuth ** 2 + 2.0 * g_thph * uuth * uuph + g_phph * uuph ** 2)
    ut = uut / alpha
    ur = uur - alpha * uut * gtr
    uth = uuth - alpha * uut * gtth
    uph = uuph - alpha * uut * gtph
    u_t = g_tt * ut + g_tr * ur + g_tth * uth + g_tph * uph
    u_r = g_tr * ut + g_rr * ur + g_rth * uth + g_rph * uph
    u_th = g_tth * ut + g_rth * ur + g_thth * uth + g_thph * uph
    u_ph = g_tph * ut + g_rph * ur + g_thph * uth + g_phph * uph
    u0 = ut
    u1 = ur / r[None,None,:]
    u2 = uth / np.pi
    u3 = uph
    u_0 = u_t
    u_1 = u_r * r[None,None,:]
    u_2 = u_th * np.pi
    u_3 = u_ph
    uu0 = alpha_alt * u0
    uu1 = u1 + alpha_alt * uu0 * g01
    uu2 = u2 + alpha_alt * uu0 * g02
    uu3 = u3 + alpha_alt * uu0 * g03

    # Calculate magnetic field
    bt = u_r * bbr + u_th * bbth + u_ph * bbph
    br = (bbr + bt * ur) / ut
    bth = (bbth + bt * uth) / ut
    bph = (bbph + bt * uph) / ut
    b_t = g_tt * bt + g_tr * br + g_tth * bth + g_tph * bph
    b_r = g_tr * bt + g_rr * br + g_rth * bth + g_rph * bph
    b_th = g_tth * bt + g_rth * br + g_thth * bth + g_thph * bph
    b_ph = g_tph * bt + g_rph * br + g_thph * bth + g_phph * bph
    b0 = bt
    b1 = br / r[None,None,:]
    b2 = bth / np.pi
    b3 = bph
    b_0 = b_t
    b_1 = b_r * r[None,None,:]
    b_2 = b_th * np.pi
    b_3 = b_ph
    bb1 = b1 * u0 - b0 * u1
    bb2 = b2 * u0 - b0 * u2
    bb3 = b3 * u0 - b0 * u3

  # Open athdf file for writing
  if kwargs['format'] == 'athdf':
    with h5py.File(kwargs['filename'], 'w') as f_out:

      # Write file-level attributes
      f_out.attrs.create('NumCycles', 0, dtype=np.int32)
      f_out.attrs.create('Time', 0.0, dtype=np.float32)
      f_out.attrs.create('Coordinates', 'kerr-schild', dtype='|S11')
      f_out.attrs.create('RootGridX1', (rf[0], rf[-1], (rf[-1] / rf[0]) ** (1.0 / len(r))),
          dtype=np.float32)
      f_out.attrs.create('RootGridX2', (thf[0], thf[-1], 1.0), dtype=np.float32)
      f_out.attrs.create('RootGridX3', (phf[0], phf[-1], 1.0), dtype=np.float32)
      f_out.attrs.create('RootGridSize', (len(r), len(th), len(ph)), dtype=np.int32)
      f_out.attrs.create('NumMeshBlocks', 1, dtype=np.int32)
      f_out.attrs.create('MeshBlockSize', (len(r), len(th), len(ph)), dtype=np.int32)
      f_out.attrs.create('MaxLevel', 0, dtype=np.int32)
      f_out.attrs.create('NumVariables', [5, 3], dtype=np.int32)
      f_out.attrs.create('DatasetNames', ['prim', 'B'], dtype='|S21')
      f_out.attrs.create('VariableNames',
          ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3'], dtype='|S21')

      # Write grid data
      f_out.create_dataset('Levels', data=(0,), dtype=np.int32)
      f_out.create_dataset('LogicalLocations', data=(0,0,0), shape=(1,3), dtype=np.int64)
      f_out.create_dataset('x1f', data=rf, shape=(1,len(rf)), dtype=np.float32)
      f_out.create_dataset('x2f', data=thf, shape=(1,len(thf)), dtype=np.float32)
      f_out.create_dataset('x3f', data=phf, shape=(1,len(phf)), dtype=np.float32)
      f_out.create_dataset('x1v', data=r, shape=(1,len(r)), dtype=np.float32)
      f_out.create_dataset('x2v', data=th, shape=(1,len(th)), dtype=np.float32)
      f_out.create_dataset('x3v', data=ph, shape=(1,len(ph)), dtype=np.float32)

      # Write cell data
      f_out.create_dataset('prim', shape=(5,1,len(ph),len(th),len(r)), dtype=np.float32)
      f_out['prim'][0,0,...] = rho
      f_out['prim'][1,0,...] = pgas
      f_out['prim'][2,0,...] = uur
      f_out['prim'][3,0,...] = uuth
      f_out['prim'][4,0,...] = uuph
      f_out.create_dataset('B', shape=(3,1,len(ph),len(th),len(r)), dtype=np.float32)
      f_out['B'][0,0,...] = bbr
      f_out['B'][1,0,...] = bbth
      f_out['B'][2,0,...] = bbph

  # Open iharm3d file for writing
  elif kwargs['format'] == 'iharm3d':
    with h5py.File(kwargs['filename'], 'w') as f_out:

      # Write header data
      f_out.create_dataset('header/version', data=('iharm-blacklight',), dtype='|S20')
      f_out.create_dataset('header/gam', data=kwargs['gamma_adi'], dtype=np.float64)
      f_out.create_dataset('header/tf', data=0.0, dtype=np.float64)
      f_out.create_dataset('header/n1', data=len(r), dtype=np.int32)
      f_out.create_dataset('header/n2', data=len(th), dtype=np.int32)
      f_out.create_dataset('header/n3', data=len(ph), dtype=np.int32)
      f_out.create_dataset('header/metric', data=('MKS',), dtype='|S20')
      f_out.create_dataset('header/n_prim', data=8, dtype=np.int32)
      f_out.create_dataset('header/prim_names', data=('RHO','UU','U1','U2','U3','B1','B2','B3'),
          dtype='|S20')
      f_out.create_dataset('header/has_electrons', data=0, dtype=np.int32)
      f_out.create_dataset('header/geom/dx1', data=dlr, dtype=np.float64)
      f_out.create_dataset('header/geom/dx2', data=dx2, dtype=np.float64)
      f_out.create_dataset('header/geom/dx3', data=dph, dtype=np.float64)
      f_out.create_dataset('header/geom/startx1', data=lrf[0], dtype=np.float64)
      f_out.create_dataset('header/geom/startx2', data=x2f[0], dtype=np.float64)
      f_out.create_dataset('header/geom/startx3', data=phf[0], dtype=np.float64)
      f_out.create_dataset('header/geom/n_dim', data=4, dtype=np.int32)
      f_out.create_dataset('header/geom/mks/r_eh', data=r_hor, dtype=np.float64)
      f_out.create_dataset('header/geom/mks/r_in', data=rf[0], dtype=np.float64)
      f_out.create_dataset('header/geom/mks/r_out', data=rf[-1], dtype=np.float64)
      f_out.create_dataset('header/geom/mks/a', data=0.0, dtype=np.float64)
      f_out.create_dataset('header/geom/mks/hslope', data=1.0, dtype=np.float64)

      # Write time
      f_out.create_dataset('t', data=0.0, dtype=np.float64)

      # Write cell data
      data = []
      data.append(rho)
      data.append(ugas)
      data.append(uu1)
      data.append(uu2)
      data.append(uu3)
      data.append(bb1)
      data.append(bb2)
      data.append(bb3)
      data = np.array(data, dtype=np.float32).transpose()
      f_out.create_dataset('prims', data=data, dtype=np.float32)

  # Open harm3d file for writing
  elif kwargs['format'] == 'harm3d':
    with open(kwargs['filename'], 'w') as f_out:

      # Write header data
      f_out.write('0.0 ')
      f_out.write('{0} {1} {2} '.format(len(r), len(th), len(ph)))
      f_out.write('{0:24.16e} {1:24.16e} {2:24.16e} '.format(lrf[0], x2f[0], phf[0]))
      f_out.write('{0:24.16e} {1:24.16e} {2:24.16e} '.format(dlr, dx2, dph))
      f_out.write('0.0 ')
      f_out.write('{0:24.16e} '.format(kwargs['gamma_adi']))
      f_out.write('{0:24.16e} '.format(rf[0]))
      f_out.write('1.0 ')
      f_out.write('8\n')

      # Write cell data
      data = []
      data.append(np.tile(lr[None,None,:], (len(ph), len(th), 1)))
      data.append(np.tile(x2[None,:,None], (len(ph), 1, len(r))))
      data.append(np.tile(ph[:,None,None], (1, len(th), len(r))))
      data.append(np.tile(r[None,None,:], (len(ph), len(th), 1)))
      data.append(np.tile(th[None,:,None], (len(ph), 1, len(r))))
      data.append(np.tile(ph[:,None,None], (1, len(th), len(r))))
      data.append(rho)
      data.append(ugas)
      data.append(u0)
      data.append(u1)
      data.append(u2)
      data.append(u3)
      data.append(b0)
      data.append(b1)
      data.append(b2)
      data.append(b3)
      data = np.array(data, dtype=np.float32).transpose()
      data.tofile(f_out)

  # Open harm3d_ext file for writing
  elif kwargs['format'] == 'harm3d_ext':
    with open(kwargs['filename'], 'w') as f_out:

      # Write header data
      f_out.write('0.0 ')
      f_out.write('{0} {1} {2} '.format(len(r), len(th), len(ph)))
      f_out.write('{0:24.16e} {1:24.16e} {2:24.16e} '.format(lrf[0], x2f[0], phf[0]))
      f_out.write('{0:24.16e} {1:24.16e} {2:24.16e} '.format(dlr, dx2, dph))
      f_out.write('0.0 ')
      f_out.write('{0:24.16e} '.format(kwargs['gamma_adi']))
      f_out.write('{0:24.16e} '.format(rf[0]))
      f_out.write('1.0 ')
      f_out.write('8\n')

      # Write cell data
      data = []
      data.append(np.tile(np.arange(len(r))[None,None,:], (len(ph), len(th), 1)))
      data.append(np.tile(np.arange(len(th))[None,:,None], (len(ph), 1, len(r))))
      data.append(np.tile(np.arange(len(ph))[:,None,None], (1, len(th), len(r))))
      data.append(np.tile(lr[None,None,:], (len(ph), len(th), 1)))
      data.append(np.tile(x2[None,:,None], (len(ph), 1, len(r))))
      data.append(np.tile(ph[:,None,None], (1, len(th), len(r))))
      data.append(np.tile(r[None,None,:], (len(ph), len(th), 1)))
      data.append(np.tile(th[None,:,None], (len(ph), 1, len(r))))
      data.append(np.tile(ph[:,None,None], (1, len(th), len(r))))
      data.append(rho)
      data.append(ugas)
      data.append(uu1)
      data.append(uu2)
      data.append(uu3)
      data.append(bb1)
      data.append(bb2)
      data.append(bb3)
      data.append(kappa)
      data.append(u0)
      data.append(u1)
      data.append(u2)
      data.append(u3)
      data.append(u_0)
      data.append(u_1)
      data.append(u_2)
      data.append(u_3)
      data.append(b0)
      data.append(b1)
      data.append(b2)
      data.append(b3)
      data.append(b_0)
      data.append(b_1)
      data.append(b_2)
      data.append(b_3)
      data.append(det_alt)
      data = np.array(data, dtype=np.float32).transpose()
      data.tofile(f_out)

  # Report invalid file format
  else:
    raise RuntimeError('Invalid format {0}.'.format(kwargs['format']))

# Parse inputs and execute main function
if __name__ == '__main__':

  # Prepare parser
  parser = argparse.ArgumentParser()

  # Prepare for filename input
  parser.add_argument('filename', help='name of simulation data file to write')
  parser.add_argument('--format', default='athdf',
     help='file format (athdf, iharm3d, harm3d, harm3d_ext) to write')

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
  parser.add_argument('--pgas_r_power', type=float, default=1.25,
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
  parser.add_argument('--Bz_R_power', type=float, default=0.625,
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

  # Prepare for miscellaneous inputs
  parser.add_argument('--gamma_adi', type=float, default=13.0/9.0,
      help='adiabatic index for formats that write internal energy')

  # Parse inputs
  args = parser.parse_args()

  # Execute main function
  main(**vars(args))
