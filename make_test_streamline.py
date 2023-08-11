#! /usr/bin/env python

"""
Script for making file with streamline points.
"""

# Python standard modules
import argparse

# Numerical modules
import numpy as np

# Main function
def main(**kwargs):

  # Parameters
  a = 0.9
  r_hor = 1.0 + (1.0 - a ** 2) ** 0.5
  center_start_r = 1.0
  center_end_r = 100.0
  th_r = lambda r: 5.0 * np.pi/180.0
  ph_r = lambda r: np.cos(2.0*np.pi * 5.0 * np.log(r / r_hor) / np.log(center_end_r / r_hor))
  num_points = 1000

  # Construct curve in spherical coordinates
  r = np.logspace(np.log10(center_start_r), np.log10(center_end_r), num_points)
  th = np.array([th_r(r_val) for r_val in r])
  ph = np.array([ph_r(r_val) for r_val in r])

  # Transform to Cartesian coordinates
  x = r * np.sin(th) * np.cos(ph)
  y = r * np.sin(th) * np.sin(ph)
  z = r * np.cos(th)

  # Write output
  with open(kwargs['filename'], 'w') as f:
    for x_val, y_val, z_val in zip(x, y, z):
      f.write('{0:24.16e} {1:24.16e} {2:24.16e}\n'.format(x_val, y_val, z_val))

# Parse inputs and execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('filename', help='name of file to write')
  args = parser.parse_args()
  main(**vars(args))
