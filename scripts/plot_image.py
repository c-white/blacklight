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
  cmap_names = ('inferno', 'PuOr', 'PiYG', 'RdBu_r')
  cmap_name_level = 'viridis'
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

  # Determine Stokes parameter
  if sum((kwargs['stokes_q'], kwargs['stokes_u'], kwargs['stokes_v'], kwargs['refinement_level'])) \
      > 1:
    raise RuntimeError('Can have at most one of Stokes Q/U/V or refinement level selected.')
  stokes = 0
  stokes = 1 if kwargs['stokes_q'] else stokes
  stokes = 2 if kwargs['stokes_u'] else stokes
  stokes = 3 if kwargs['stokes_v'] else stokes
  stokes = None if kwargs['refinement_level'] else stokes

  # Prepare metadata
  width_rg = kwargs['width']
  mass_msun = kwargs['mass']
  distance_pc = kwargs['distance']
  max_level = kwargs['max_level']

  # Read data from .npz file
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:

      # Read root image
      polarization = f['image'].shape[0] == 4
      if not polarization and stokes != 0:
        raise RuntimeError('No polarization data in file.')
      if stokes is not None:
        image = f['image'][stokes,:,:]
      else:
        image = np.zeros_like(f['image'][0,:,:])

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
          if stokes is not None:
            image_adaptive[level] = f['adaptive_image_{0}'.format(level)][stokes,:,:,:]
          else:
            image_adaptive[level] = \
                level * np.ones_like(f['adaptive_image_{0}'.format(level)][0,:,:,:])

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
    if stokes is None:
      raise RuntimeError('Adaptive refinement not supported by npy format.')
    image = np.load(kwargs['filename_data'])
    polarization = image.shape[0] == 4
    if not polarization and stokes != 0:
      raise RuntimeError('No polarization data in file.')
    image = image[stokes,:,:]
    max_level = 0

  # Read image data from raw file
  else:
    if stokes is None:
      raise RuntimeError('Adaptive refinement not supported by raw format.')
    if stokes != 0:
      raise RuntimeError('No polarization data in file.')
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

  # Calculate intensity scale
  if stokes is not None:
    vmax = np.nanmax(np.abs(image))
    for level in range(1, max_level + 1):
      vmax_adaptive = np.nanmax(np.abs(image_adaptive[level]))
      vmax = np.nanmax((vmax, vmax_adaptive))
    if stokes == 0:
      vmin = 0.0
    else:
      vmin = -vmax
    scale_exponent = int('{0:24.16e}'.format(vmax).split('e')[1])
    scale = 10.0 ** scale_exponent
    vmin /= scale
    vmax /= scale
    image /= scale
    for level in range(1, max_level + 1):
      image_adaptive[level] /= scale
    stokes_str = 'I'
    stokes_str = 'Q' if kwargs['stokes_q'] else stokes_str
    stokes_str = 'U' if kwargs['stokes_u'] else stokes_str
    stokes_str = 'V' if kwargs['stokes_v'] else stokes_str
    if scale_exponent == 0:
      label = r'$' + stokes_str + r'_\nu$ ($\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
    elif scale_exponent == 1:
      label = r'$' + stokes_str + r'_\nu$ ($10\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
    else:
      label = r'$' + stokes_str + r'_\nu$ ($10^{' + repr(scale_exponent) \
          + r'}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
    tick_locs = None

  # Calculate refinement scale
  else:
    vmin = -0.5
    vmax = max_level + 0.5
    label = 'level'
    tick_locs = np.arange(0, max_level + 1)

  # Define colormap
  if stokes is not None:
    cmap = plt.get_cmap(cmap_names[stokes])
    cmap.set_bad(nan_color)
    bounds = None
  else:
    cmap_continuous = plt.get_cmap(cmap_name_level)
    cmap_vals = [cmap_continuous(float(x) / max_level) for x in tick_locs]
    cmap = LinearSegmentedColormap.from_list('levels', cmap_vals, max_level + 1)
    bounds = np.linspace(vmin, vmax, max_level + 2)

  # Plot root image
  plt.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal', origin='lower', extent=extent,
      interpolation=interpolation)

  # Plot adaptive image
  if max_level > 0:
    for level in range(1, max_level + 1):
      for block in range(num_blocks[level]):
        plt.imshow(image_adaptive[level][block,:,:], cmap=cmap, vmin=vmin, vmax=vmax,
            aspect='equal', origin='lower', extent=extent_adaptive[level][block,:],
            interpolation=interpolation)

  # Make colorbar
  cb = plt.colorbar(ticks=tick_locs, boundaries=bounds)
  cb.set_label(label, labelpad=labelpad)
  if stokes is None:
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
  parser.add_argument('-r', '--refinement_level', action='store_true',
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
