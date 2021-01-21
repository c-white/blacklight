#! /usr/bin/env python

"""
Script for calculating total flux from outputs produced by blacklight.
"""

# Python standard modules
import argparse

# Other modules
import numpy as np

# Main function
def main(**kwargs):

  # Parameters
  c = 2.99792458e10
  gg_msun = 1.32712440018e26
  kb = 1.380649e-16
  pc = 9.69394202136e18 / np.pi
  jy = 1.0e-23
  data_format = np.float64

  # Read image data from NumPy file
  if kwargs['image_data_file'][-4:] == '.npy':
    tb = np.load(kwargs['image_data_file'])

  # Read image data from raw file
  else:
    image = np.fromfile(kwargs['image_data_file'], dtype=data_format)

  # Calculate flux
  d = kwargs['distance'] * pc
  w = kwargs['width'] * gg_msun * kwargs['mass'] / c**2
  inu = 2.0 * kwargs['frequency']**2 * kb / c**2 * tb
  fnu = np.mean(inu) * (w / d) ** 2
  fnu_jy = fnu / jy

  # Report results
  print('F_nu = {0} Jy'.format(repr(fnu_jy)))

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('image_data_file', help='file containing raw image data')
  parser.add_argument('mass', type=float, help='black hole mass in solar masses')
  parser.add_argument('distance', type=float, help='distance to black hole in parsecs')
  parser.add_argument('width', type=float, help='full width of figure in gravitational radii')
  parser.add_argument('frequency', type=float, help='observed frequency in Hz')
  args = parser.parse_args()
  main(**vars(args))
