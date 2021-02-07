#! /usr/bin/env python

"""
Script for plotting outputs produced by blacklight.
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
  data_format = np.float64
  cmap_name = 'inferno'
  nan_color = 'gray'
  dpi = 300

  # Define colormap
  cmap = plt.get_cmap(cmap_name)
  cmap.set_bad(nan_color)

  # Read image data from NumPy file
  if kwargs['image_data_file'][-4:] == '.npy':
    image = np.load(kwargs['image_data_file'])
  elif kwargs['image_data_file'][-4:] == '.npz':
    image = np.load(kwargs['image_data_file'])['image']

  # Read image data from raw file
  else:
    image = np.fromfile(kwargs['image_data_file'], dtype=data_format)
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
  plt.imshow(image, origin='lower', cmap=cmap, vmin=0.0)
  plt.colorbar()
  plt.tight_layout()
  plt.savefig(kwargs['image_plot_file'], dpi=dpi)

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('image_data_file', help='file containing raw image data')
  parser.add_argument('image_plot_file', help='image file to create')
  args = parser.parse_args()
  main(**vars(args))
