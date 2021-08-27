#! /usr/bin/env python

"""
Script for calculating total flux from outputs produced by blacklight.
"""

# Python standard modules
import argparse

# Numerical modules
import numpy as np

# Main function
def main(**kwargs):

  # Parameters
  c = 2.99792458e10
  gg_msun = 1.32712440018e26
  pc = 9.69394202136e18 / np.pi
  jy = 1.0e-23
  data_format = np.float64

  # Prepare metadata
  width_rg = kwargs['width']
  mass_msun = kwargs['mass']
  distance_pc = kwargs['distance']
  max_level = kwargs['max_level']

  # Read data from .npz file
  if kwargs['filename_data'][-4:] == '.npz':
    with np.load(kwargs['filename_data']) as f:

      # Read root image
      image = f['image'][:]
      polarization = image.shape[0] == 4

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
          image_adaptive[level] = f['adaptive_image_{0}'.format(level)][:]

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
    image = np.load(kwargs['filename_data'])
    polarization = image.shape[0] == 4
    max_level = 0

  # Read image data from raw file
  else:
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
    polarization = False
    max_level = 0

  # Calculate image size
  if width_rg is None:
    raise RuntimeError('Must supply width.')
  if mass_msun is None:
    raise RuntimeError('Must supply mass.')
  if distance_pc is None:
    raise RuntimeError('Must supply distance.')
  rg = gg_msun * mass_msun / c ** 2
  width = 2.0 * np.arctan(0.5 * width_rg * rg / (distance_pc * pc))

  # Calculate flux without adaptive refinement
  if max_level == 0:
    flux = np.mean(image, axis=(1,2)) * width ** 2 / jy

  # Calculate flux with adaptive refinement
  else:

    # Disassemble root image into blocks
    block_res = image_adaptive[1].shape[-1]
    num_blocks_root_linear = image.shape[-1] / block_res
    image_adaptive[0] = np.reshape(image,
        (-1, num_blocks_root_linear, block_res, num_blocks_root_linear, block_res))
    image_adaptive[0] = np.swapaxes(image_adaptive[0], 2, 3)
    image_adaptive[0] = np.reshape(image_adaptive[0], (-1, num_blocks[0], block_res, block_res))

    # Describe disassembled image locations
    block_locs[0] = np.empty((num_blocks[0], 2), dtype=int)
    for block in range(num_blocks[0]):
      block_locs[0][block,1] = block % num_blocks_root_linear
      block_locs[0][block,0] = block / num_blocks_root_linear

    # Prepare flags indicating blocks to be counted
    flags = {}
    for level in range(max_level + 1):
      flags[level] = [True] * num_blocks[level]

    # Determine which blocks are masked by those at higher levels
    for level in range(1, max_level + 1):
      for block_fine in range(num_blocks[level]):
        x_loc_fine = block_locs[level][block_fine,1]
        y_loc_fine = block_locs[level][block_fine,0]
        x_loc_coarse = x_loc_fine / 2
        y_loc_coarse = y_loc_fine / 2
        for block_coarse in range(num_blocks[level-1]):
          x_loc_coarse_test = block_locs[level-1][block_coarse,1]
          y_loc_coarse_test = block_locs[level-1][block_coarse,0]
          if x_loc_coarse == x_loc_coarse_test and y_loc_coarse == y_loc_coarse_test:
            flags[level-1][block_coarse] = False

    # Prepare fluxes
    if polarization:
      flux = np.zeros(4)
    else:
      flux = np.zeros(1)

    # Calculate flux from non-masked blocks
    block_width_root = width / num_blocks_root_linear
    for level in range(max_level + 1):
      block_width = block_width_root / 2 ** level
      for block in range(num_blocks[level]):
        if flags[level][block]:
          flux += np.mean(image_adaptive[level][:,block,:,:], axis=(1,2)) * block_width ** 2
    flux /= jy

  # Report results
  print('')
  if polarization:
    print('I: F_nu = {0} Jy'.format(repr(flux[0])))
    print('Q: F_nu = {0} Jy'.format(repr(flux[1])))
    print('U: F_nu = {0} Jy'.format(repr(flux[2])))
    print('V: F_nu = {0} Jy'.format(repr(flux[3])))
  else:
    print('F_nu = {0} Jy'.format(repr(flux[0])))
  print('')

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('filename_data', help='name of file containing raw image data')
  parser.add_argument('-w', '--width', type=float,
      help='full width of image in gravitational radii')
  parser.add_argument('-m', '--mass', type=float, help='black hole mass in solar masses')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in parsecs')
  parser.add_argument('-l', '--max_level', type=int,
      help='maximum adaptive level to use in calculation')
  args = parser.parse_args()
  main(**vars(args))
