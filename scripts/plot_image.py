#! /usr/bin/env python

"""
Script for plotting outputs produced by ray_trace.
"""

# Python standard modules
import argparse

# Python plotting modules
import matplotlib
matplotlib.use('agg')
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt

# Other modules
import numpy as np

# Main function
def main(**kwargs):

  # Parameters
  data_format = np.float32
  dpi = 300

  # Read image data
  image = np.fromfile(kwargs['image_data_file'], dtype=data_format)

  # Reshape image
  pix_num = len(image)
  pix_res = int(pix_num ** 0.5)
  if pix_res ** 2 == pix_num:
    pass
  elif (pix_res + 1) ** 2 == pix_num:
    pix_res += 1
  else:
    raise RuntimeError('Image data not square')
  image = np.reshape(image, (pix_res, pix_res))

  # Plot figure
  plt.imshow(image, origin='lower')
  plt.colorbar()
  plt.savefig(kwargs['image_plot_file'], dpi=dpi)

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('image_data_file', help='file containing raw image data')
  parser.add_argument('image_plot_file', help='image file to create')
  args = parser.parse_args()
  main(**vars(args))