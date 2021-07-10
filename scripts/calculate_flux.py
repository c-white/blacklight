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
  pc = 9.69394202136e18 / np.pi
  jy = 1.0e-23
  data_format = np.float64

  # Prepare to read data
  mass = None
  distance = None
  width = None

  # Read image data from .npz file
  if kwargs['image_data_file'][-4:] == '.npz':
    data = np.load(kwargs['image_data_file'])
    image = data['image']
    if 'image_width' in data.keys():
      width = data['image_width'][0]
    if 'mass_msun' in data.keys():
      mass = data['mass_msun'][0]
    if len(image.shape) == 3:
      image = image[0,:,:]

  # Read image data from .npy file
  elif kwargs['image_data_file'][-4:] == '.npy':
    image = np.load(kwargs['image_data_file'])
    if len(image.shape) == 3:
      image = image[0,:,:]

  # Read image data from raw file
  else:
    image = np.fromfile(kwargs['image_data_file'], dtype=data_format)

  # Check for mass
  if kwargs['mass'] is None and mass is None:
    raise RuntimeError('Must supply mass.')
  elif kwargs['mass'] is None and mass is not None:
    pass
  elif kwargs['mass'] is not None and mass is None:
    mass = kwargs['mass']
  else:
    if not np.isclose(kwargs['mass'], mass, rtol=1.0e-5, atol=0.0):
      print('Warning: Using "mass" that does not match file "mass_msun".')
    mass = kwargs['mass']

  # Check for distance
  if kwargs['distance'] is None and distance is None:
    raise RuntimeError('Must supply distance.')
  elif kwargs['distance'] is None and distance is not None:
    pass
  elif kwargs['distance'] is not None and distance is None:
    distance = kwargs['distance']
  else:
    if not np.isclose(kwargs['distance'], distance, rtol=1.0e-5, atol=0.0):
      print('Warning: Using "distance" that does not match file.')
    distance = kwargs['distance']

  # Check for width
  if kwargs['width'] is None and width is None:
    raise RuntimeError('Must supply width.')
  elif kwargs['width'] is None and width is not None:
    pass
  elif kwargs['width'] is not None and width is None:
    width = kwargs['width']
  else:
    if not np.isclose(kwargs['width'], width, rtol=1.0e-5, atol=0.0):
      print('Warning: Using "width" that does not match file "image_width".')
    width = kwargs['width']

  # Calculate flux
  d = distance * pc
  w = width * gg_msun * mass / c**2
  fnu = np.mean(image) * (w / d) ** 2
  fnu_jy = fnu / jy

  # Report results
  print('F_nu = {0} Jy'.format(repr(fnu_jy)))

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('image_data_file', help='file containing raw image data')
  parser.add_argument('-m', '--mass', type=float, help='black hole mass in solar masses')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in parsecs')
  parser.add_argument('-w', '--width', type=float,
      help='full width of figure in gravitational radii')
  args = parser.parse_args()
  main(**vars(args))
