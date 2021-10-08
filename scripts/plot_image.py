#! /usr/bin/env python

"""
Script for plotting outputs produced by blacklight.
"""

# Python standard modules
import argparse

# Numerical modules
import numpy as np

# Main function
def main(**kwargs):

  # Parameters
  c = 2.99792458e10;
  gg_msun = 1.32712440018e26;
  pc = 9.69394202136e18 / np.pi
  muas = np.pi / 180.0 / 60.0 / 60.0 / 1.0e6
  data_format = np.float64

  # Plotting parameters
  cmap_i = 'inferno'
  cmap_q = 'PuOr'
  cmap_u = 'PiYG'
  cmap_v = 'RdBu_r'
  cmap_level = 'plasma'
  cmap_name = 'viridis'
  nan_color = 'gray'
  interpolation = 'none'
  labelpad = 5
  dpi = 300

  # Import plotting module
  import matplotlib
  matplotlib.use('agg')
  if not kwargs['notex']:
    matplotlib.rc('text', usetex=True)
  import matplotlib.pyplot as plt
  from matplotlib.colors import LinearSegmentedColormap

  # Verify input
  if sum((kwargs['stokes_q'], kwargs['stokes_u'], kwargs['stokes_v'], kwargs['name'] is not None,
      kwargs['rendering'] is not None, kwargs['refinement_level'])) > 1:
    raise RuntimeError('Can have at most one of Stokes Q/U/V, named quantity, rendering, or' \
        + ' refinement level selected.')

  # Prepare metadata
  width_rg = kwargs['width']
  mass_msun = kwargs['mass']
  distance_pc = kwargs['distance']
  max_level = kwargs['max_level']

  # Read data from .npz file
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:

      # Read root image
      if kwargs['refinement_level']:
        try:
          image = np.zeros_like(f['I_nu'][:])
        except KeyError:
          raise RuntimeError('No intensity data in file, so image cannot be refined.')
      elif kwargs['rendering'] is not None:
        try:
          image = f['rendering'][kwargs['rendering']-1,...]
        except:
          raise RuntimeError('Rendering {0} not found in file.'.format(kwargs['rendering']))
      elif kwargs['name'] is not None:
        try:
          image = f[kwargs['name']][:]
        except:
          raise RuntimeError('{0} not found in file.'.format(kwargs['name']))
      elif kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']:
        try:
          if kwargs['stokes_q']:
            image = f['Q_nu'][:]
          elif kwargs['stokes_u']:
            image = f['U_nu'][:]
          else:
            image = f['V_nu'][:]
        except KeyError:
          raise RuntimeError('No polarization data in file.')
      else:
        try:
          image = f['I_nu'][:]
        except KeyError:
          raise RuntimeError('No intensity data in file.')

      # Read adaptive image
      if max_level is None:
        max_level = f['adaptive_num_levels'][0]
      else:
        max_level = min(max_level, f['adaptive_num_levels'][0])
      if max_level > 0:
        num_blocks = {}
        block_locs = {}
        image_adaptive = {}
        num_blocks[0] = f['adaptive_num_blocks'][0]
        for level in range(1, max_level + 1):
          num_blocks[level] = f['adaptive_num_blocks'][level]
          block_locs[level] = f['adaptive_block_locs_{0}'.format(level)][:]
          if kwargs['refinement_level']:
            image_adaptive[level] = \
                level * np.ones_like(f['adaptive_I_nu_{0}'.format(level)])
          elif kwargs['rendering'] is not None:
            image_adaptive[level] = \
                f['adaptive_rendering_{0}'.format(level)][kwargs['rendering']-1,...]
          elif kwargs['name'] is not None:
            image_adaptive[level] = f['adaptive_{0}_{1}'.format(kwargs['name'], level)][:]
          elif kwargs['stokes_q']:
            image_adaptive[level] = f['adaptive_Q_nu_{0}'.format(level)][:]
          elif kwargs['stokes_u']:
            image_adaptive[level] = f['adaptive_U_nu_{0}'.format(level)][:]
          elif kwargs['stokes_v']:
            image_adaptive[level] = f['adaptive_V_nu_{0}'.format(level)][:]
          else:
            image_adaptive[level] = f['adaptive_I_nu_{0}'.format(level)][:]

      # Read metadata
      if 'width' in f.keys():
        if width_rg is None:
          width_rg = f['width'][0]
        elif not np.isclose(f['width'][0], width_rg):
          raise RuntimeError('Input width {0} does not match file value {1}.'.format(width_rg,
              f['width'][0]))
      if 'mass_msun' in f.keys():
        if mass_msun is None:
          mass_msun = f['mass_msun'][0]
        elif not np.isclose(f['mass_msun'][0], mass_msun):
          raise RuntimeError('Input mass {0} does not match file value {1}.'.format(mass_msun,
              f['mass_msun'][0]))

  # Read image data from .npy file
  elif kwargs['filename_data'][-4:] == '.npy':
    if kwargs['refinement_level']:
      raise RuntimeError('Adaptive refinement not supported by npy format.')
    if kwargs['rendering'] is not None:
      raise RuntimeError('Renderings not supported by npy format.')
    if kwargs['name'] is not None:
      raise RuntimeError('Arbitrary named quantities not supported by npy format.')
    image = np.load(kwargs['filename_data'])
    polarization = image.shape[0] == 4
    if not polarization and (kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']):
      raise RuntimeError('No polarization data in file.')
    if kwargs['stokes_q']:
      image = image[1,...]
    elif kwargs['stokes_u']:
      image = image[2,...]
    elif kwargs['stokes_v']:
      image = image[3,...]
    else:
      image = image[0,...]
    max_level = 0

  # Read image data from raw file
  else:
    if kwargs['refinement_level']:
      raise RuntimeError('Adaptive refinement not supported by raw format.')
    if kwargs['rendering'] is not None:
      raise RuntimeError('Renderings not supported by raw format.')
    if kwargs['name'] is not None:
      raise RuntimeError('Arbitrary named quantities not supported by raw format.')
    if kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']:
      raise RuntimeError('Polarization not supported by raw format.')
    image = np.fromfile(kwargs['filename_data'], dtype=data_format)
    pix_num = len(image)
    pix_res = int(pix_num ** 0.5)
    if pix_res ** 2 == pix_num:
      pass
    elif (pix_res + 1) ** 2 == pix_num:
      pix_res += 1
    else:
      raise RuntimeError('Image data not square.')
    image = np.reshape(image, (pix_res, pix_res))
    max_level = 0

  # Calculate root grid in pixels
  if kwargs['axes'] == 'pixel' or (kwargs['axes'] is None and width_rg is None):
    extent = np.array((-0.5, image.shape[0] - 0.5, -0.5, image.shape[1] - 0.5))
    x_label = r'$x$-pixel'
    y_label = r'$y$-pixel'

  # Calculate root grid in gravitational radii
  elif kwargs['axes'] == 'rg' or (kwargs['axes'] is None and mass_msun is None):
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    half_width = 0.5 * width_rg
    scale_exponent = int('{0:24.16e}'.format(half_width).split('e')[1])
    if scale_exponent in (0, 1):
      scale = 1.0
      x_label = r'$x$ ($GM/c^2$)'
      y_label = r'$y$ ($GM/c^2$)'
    else:
      scale = 10.0 ** scale_exponent
      x_label = r'$x$ ($10^{' + repr(scale_exponent) + r'}\ GM/c^2$)'
      y_label = r'$y$ ($10^{' + repr(scale_exponent) + r'}\ GM/c^2$)'
    half_width /= scale
    extent = np.array((-half_width, half_width, -half_width, half_width))

  # Calculate root grid in cm
  elif kwargs['axes'] == 'cm' or (kwargs['axes'] is None and distance_pc is None):
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    if mass_msun is None:
      raise RuntimeError('Must supply mass.')
    rg = gg_msun * mass_msun / c ** 2
    half_width = 0.5 * width_rg * rg
    scale_exponent = int('{0:24.16e}'.format(half_width).split('e')[1])
    if scale_exponent in (0, 1):
      scale = 1.0
      x_label = r'$x$ ($\mathrm{cm}$)'
      y_label = r'$y$ ($\mathrm{cm}$)'
    else:
      scale = 10.0 ** scale_exponent
      x_label = r'$x$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{cm}$)'
      y_label = r'$y$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{cm}$)'
    half_width /= scale
    extent = np.array((-half_width, half_width, -half_width, half_width))

  # Calculate root grid in muas
  else:
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    if mass_msun is None:
      raise RuntimeError('Must supply mass.')
    if distance_pc is None:
      raise RuntimeError('Must supply distance.')
    rg = gg_msun * mass_msun / c ** 2
    half_width = np.arctan(0.5 * width_rg * rg / (distance_pc * pc)) / muas
    scale_exponent = int('{0:24.16e}'.format(half_width).split('e')[1])
    if scale_exponent in (0, 1):
      scale = 1.0
      x_label = r'$x$ ($\mathrm{\mu as}$)'
      y_label = r'$y$ ($\mathrm{\mu as}$)'
    else:
      scale = 10.0 ** scale_exponent
      x_label = r'$x$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{\mu as}$)'
      y_label = r'$y$ ($10^{' + repr(scale_exponent) + r'}\ \mathrm{\mu as}$)'
    half_width /= scale
    extent = np.array((-half_width, half_width, -half_width, half_width))

  # Calculate adaptive grid
  if max_level > 0:
    num_blocks_root_linear = image.shape[-1] / image_adaptive[1].shape[-1]
    block_width = (extent[1] - extent[0]) / num_blocks_root_linear
    extent_adaptive = {}
    for level in range(1, max_level + 1):
      block_width_level = block_width / 2 ** level
      extent_adaptive[level] = np.empty((num_blocks[level], 4))
      for block in range(num_blocks[level]):
        x_loc = block_locs[level][block,1]
        y_loc = block_locs[level][block,0]
        extent_adaptive[level][block,0] = extent[0] + x_loc * block_width_level
        extent_adaptive[level][block,1] = extent[0] + (x_loc + 1) * block_width_level
        extent_adaptive[level][block,2] = extent[2] + y_loc * block_width_level
        extent_adaptive[level][block,3] = extent[2] + (y_loc + 1) * block_width_level

  # Calculate refinement parameters
  if kwargs['refinement_level']:
    label = 'level'
    vmin = -0.5
    vmax = max_level + 0.5
    tick_locs = np.arange(0, max_level + 1)
    cmap_continuous = plt.get_cmap(cmap_level)
    cmap_vals = [cmap_continuous(float(x) / max_level) for x in tick_locs]
    cmap = LinearSegmentedColormap.from_list('levels', cmap_vals, max_level + 1)
    bounds = np.linspace(vmin, vmax, max_level + 2)

  # Calculate rendering parameters
  elif kwargs['rendering'] is not None:
    r_vals = 3.2406 * image[0,...] - 1.5372 * image[1,...] - 0.4986 * image[2,...]
    g_vals = -0.9689 * image[0,...] + 1.8758 * image[1,...] + 0.0415 * image[2,...]
    b_vals = 0.0557 * image[0,...] - 0.204 * image[1,...] + 1.057 * image[2,...]
    r_vals = np.clip(r_vals, 0.0, 1.0)
    g_vals = np.clip(g_vals, 0.0, 1.0)
    b_vals = np.clip(b_vals, 0.0, 1.0)
    r_vals = np.where(r_vals <= 0.0031308, 12.92 * r_vals, 1.055 * r_vals ** (1.0 / 2.4) - 0.055)
    g_vals = np.where(g_vals <= 0.0031308, 12.92 * g_vals, 1.055 * g_vals ** (1.0 / 2.4) - 0.055)
    b_vals = np.where(b_vals <= 0.0031308, 12.92 * b_vals, 1.055 * b_vals ** (1.0 / 2.4) - 0.055)
    image = np.concatenate((r_vals[:,:,None], g_vals[:,:,None], b_vals[:,:,None]), axis=2)
    for level in range(1, max_level + 1):
      r_vals = 3.2406 * image_adpative[level][0,...] - 1.5372 * image_adaptive[level][1,...] \
          - 0.4986 * image_adaptive[level][2,...]
      g_vals = -0.9689 * image_adaptive[level][0,...] + 1.8758 * image_adaptive[level][1,...] \
          + 0.0415 * image_adaptive[level][2,...]
      b_vals = 0.0557 * image_adaptive[level][0,...] - 0.204 * image_adaptive[level][1,...] \
          + 1.057 * image_adaptive[level][2,...]
      r_vals = np.clip(r_vals, 0.0, 1.0)
      g_vals = np.clip(g_vals, 0.0, 1.0)
      b_vals = np.clip(b_vals, 0.0, 1.0)
      r_vals = np.where(r_vals <= 0.0031308, 12.92 * r_vals, 1.055 * r_vals ** (1.0 / 2.4) - 0.055)
      g_vals = np.where(g_vals <= 0.0031308, 12.92 * g_vals, 1.055 * g_vals ** (1.0 / 2.4) - 0.055)
      b_vals = np.where(b_vals <= 0.0031308, 12.92 * b_vals, 1.055 * b_vals ** (1.0 / 2.4) - 0.055)
      image_adaptive[level] = \
          np.concatenate((r_vals[:,:,:,None], g_vals[:,:,:,None], b_vals[:,:,:,None]), axis=3)

  # Calculate named quantity parameters
  elif kwargs['name'] is not None:
    label_mid = {'time': r'\mathrm{min}(t)', 'length': r'l', 'lambda': r'\lambda',
        'emission': r'\int j_\nu \nu^{-2}\ \mathrm{d}\lambda', 'tau': r'\tau', 'rho': r'\rho',
        'n_e': r'n_\mathrm{e}', 'p_gas': r'p_\mathrm{gas}', 'Theta_e': r'\Theta_\mathrm{e}',
        'B': r'B', 'sigma': r'\sigma', 'beta_inverse': r'\beta^{-1}'}
    label_unit = {'time': r' ($\mathrm{s}$)', 'length': r' ($\mathrm{cm}$)',
        'lambda': r' ($\mathrm{s\ cm}$)',
        'emission': r' ($\mathrm{erg\ s^2\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)', 'tau': r'',
        'rho': r' ($\mathrm{g\ cm^{-3}}$)', 'n_e': r' ($\mathrm{cm^{-3}}$)',
        'p_gas': r' ($\mathrm{dyne\ cm^{-2}}$)', 'Theta_e': r'', 'B': r' ($\mathrm{gauss}$)',
        'sigma': r'', 'beta_inverse': r''}
    try:
      if kwargs['name'] in ('time', 'length', 'lambda', 'tau', 'emission'):
        key_mid = kwargs['name']
        label = r'$' + label_mid[kwargs['name']] + r'$'
      elif kwargs['name'][:11] == 'lambda_ave_':
        key_mid = kwargs['name'][11:]
        label = r'$\langle ' + label_mid[key_mid] + r' \rangle_\lambda$'
      elif kwargs['name'][:13] == 'emission_ave_':
        key_mid = kwargs['name'][13:]
        label = r'$\langle ' + label_mid[key_mid] + r' \rangle_{j_\nu / \nu^2}$'
      elif kwargs['name'][:8] == 'tau_int_':
        key_mid = kwargs['name'][8:]
        label = r'$\langle ' + label_mid[key_mid] + r' \rangle_\tau$'
      else:
        raise KeyError
      label += label_unit[key_mid]
    except KeyError:
      underscore = '_' if kwargs['notex'] else r'\_'
      label = ''.join([c if c != '_' else underscore for c in kwargs['name']])
    vmin = np.nanmin(image)
    vmax = np.nanmax(image)
    for level in range(1, max_level + 1):
      vmin_adaptive = np.nanmin(image_adaptive[level])
      vmax_adaptive = np.nanmax(image_adaptive[level])
      vmin = np.nanmin((vmin, vmin_adaptive))
      vmax = np.nanmax((vmax, vmax_adaptive))
    tick_locs = None
    cmap = plt.get_cmap(cmap_name)
    cmap.set_bad(nan_color)
    bounds = None

  # Calculate intensity parameters
  else:
    vmax = np.nanmax(np.abs(image))
    for level in range(1, max_level + 1):
      vmax_adaptive = np.nanmax(np.abs(image_adaptive[level]))
      vmax = np.nanmax((vmax, vmax_adaptive))
    if kwargs['stokes_q'] or kwargs['stokes_u'] or kwargs['stokes_v']:
      vmin = -vmax
    else:
      vmin = 0.0
    scale_exponent = int('{0:24.16e}'.format(vmax).split('e')[1])
    scale = 10.0 ** scale_exponent
    vmin /= scale
    vmax /= scale
    image /= scale
    for level in range(1, max_level + 1):
      image_adaptive[level] /= scale
    if kwargs['stokes_q']:
      stokes_str = 'Q'
      cmap = plt.get_cmap(cmap_q)
    elif kwargs['stokes_u']:
      stokes_str = 'U'
      cmap = plt.get_cmap(cmap_u)
    elif kwargs['stokes_v']:
      stokes_str = 'V'
      cmap = plt.get_cmap(cmap_v)
    else:
      stokes_str = 'I'
      cmap = plt.get_cmap(cmap_i)
    if scale_exponent == 0:
      label = r'$' + stokes_str + r'_\nu$ ($\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
    elif scale_exponent == 1:
      label = r'$' + stokes_str + r'_\nu$ ($10\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
    else:
      label = r'$' + stokes_str + r'_\nu$ ($10^{' + repr(scale_exponent) \
          + r'}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
    tick_locs = None
    cmap.set_bad(nan_color)
    bounds = None

  # Plot root image
  if kwargs['rendering'] is not None:
    plt.imshow(image, aspect='equal', origin='lower', extent=extent, interpolation=interpolation)
  else:
    plt.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal', origin='lower',
        extent=extent, interpolation=interpolation)

  # Plot adaptive image
  for level in range(1, max_level + 1):
    for block in range(num_blocks[level]):
      if kwargs['rendering'] is not None:
        plt.imshow(image_adaptive[level][block,...], aspect='equal', origin='lower',
            extent=extent_adaptive[level][block,:], interpolation=interpolation)
      else:
        plt.imshow(image_adaptive[level][block,...], cmap=cmap, vmin=vmin, vmax=vmax,
            aspect='equal', origin='lower', extent=extent_adaptive[level][block,:],
            interpolation=interpolation)

  # Make colorbar
  if kwargs['rendering'] is None:
    cb = plt.colorbar(ticks=tick_locs, boundaries=bounds)
    cb.set_label(label, labelpad=labelpad)
    if kwargs['refinement_level']:
      cb.ax.tick_params(length=0)

  # Adjust axes
  plt.xlim(extent[0], extent[1])
  plt.xlabel(x_label)
  plt.ylim(extent[2], extent[3])
  plt.ylabel(y_label)

  # Adjust layout
  plt.tight_layout()

  # Save figure
  plt.savefig(kwargs['filename_plot'], dpi=dpi)

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('filename_data', help='name of file containing raw image data')
  parser.add_argument('filename_plot', help='name of image file to write')
  parser.add_argument('-q', '--stokes_q', action='store_true',
      help='flag indicating Stokes Q_nu should be plotted')
  parser.add_argument('-u', '--stokes_u', action='store_true',
      help='flag indicating Stokes U_nu should be plotted')
  parser.add_argument('-v', '--stokes_v', action='store_true',
      help='flag indicating Stokes V_nu should be plotted')
  parser.add_argument('-n', '--name', help='name of quantity to be plotted')
  parser.add_argument('-r', '--rendering', type=int, help='number of rendering to be plotted')
  parser.add_argument('--refinement_level', action='store_true',
      help='flag indicating refinement level should be plotted')
  parser.add_argument('-w', '--width', type=float,
      help='full width of figure in gravitational radii')
  parser.add_argument('-m', '--mass', type=float, help='black hole mass in solar masses')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in parsecs')
  parser.add_argument('-a', '--axes', choices=('pixel','rg','cm','muas'), help='axes labels')
  parser.add_argument('-l', '--max_level', type=int, help='maximum adaptive level to plot')
  parser.add_argument('--notex', action='store_true',
      help='flag indicating external tex distribution should not be used')
  args = parser.parse_args()
  main(**vars(args))
