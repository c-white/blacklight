#! /usr/bin/env python

"""
Script for plotting rendering outputs produced by Blacklight.
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

  # Plotting parameters
  interpolation = 'none'
  dpi = 300

  # Import plotting module
  import matplotlib
  matplotlib.use('agg')
  if not kwargs['notex']:
    matplotlib.rc('text', usetex=True)
  import matplotlib.pyplot as plt

  # Prepare metadata
  distance_pc = kwargs['distance']
  max_level = kwargs['max_level']

  # Read data from .npz file
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:

      # Read metadata
      width_rg = f['width'][0]
      mass_msun = f['mass_msun'][0]

      # Read root image
      try:
        image = f['rendering'][kwargs['rendering']-1,...]
      except:
        raise RuntimeError('Rendering {0} not found in file.'.format(kwargs['rendering']))

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
          image_adaptive[level] = \
              f['adaptive_rendering_{0}'.format(level)][kwargs['rendering']-1,...]

  # Read image data from .npy file
  elif kwargs['filename_data'][-4:] == '.npy':
    raise RuntimeError('Renderings not supported by npy format.')

  # Read image data from raw file
  else:
    raise RuntimeError('Renderings not supported by raw format.')

  # Calculate root grid in pixels
  if kwargs['axes'] == 'pixel':
    extent = np.array((-0.5, image.shape[1] - 0.5, -0.5, image.shape[2] - 0.5))
    x_label = r'$x$-pixel'
    y_label = r'$y$-pixel'

  # Calculate root grid in gravitational radii
  elif kwargs['axes'] == 'rg':
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

  # Calculate sRGB1 colors from XYZ1 colors
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

  # Plot root image
  plt.imshow(image, aspect='equal', origin='lower', extent=extent, interpolation=interpolation)

  # Plot adaptive image
  for level in range(1, max_level + 1):
    for block in range(num_blocks[level]):
      plt.imshow(image_adaptive[level][block,...], aspect='equal', origin='lower',
          extent=extent_adaptive[level][block,:], interpolation=interpolation)

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
  parser.add_argument('rendering', type=int, help='number of rendering to be plotted')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in parsecs')
  parser.add_argument('-a', '--axes', choices=('pixel','rg','cm','muas'), help='axes labels')
  parser.add_argument('-l', '--max_level', type=int, help='maximum adaptive level to plot')
  parser.add_argument('--notex', action='store_true',
      help='flag indicating external tex distribution should not be used')
  args = parser.parse_args()
  main(**vars(args))
