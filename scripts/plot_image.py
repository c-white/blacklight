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
  c = 2.99792458e10;
  gg_msun = 1.32712440018e26;
  pc = 9.69394202136e18 / np.pi
  data_format = np.float64
  cmap_name = 'inferno'
  nan_color = 'gray'
  interpolation = 'none'
  dpi = 300

  # Define colormap
  cmap = plt.get_cmap(cmap_name)
  cmap.set_bad(nan_color)

  # Prepare to read data
  im_camera = None
  im_width = None
  mass_msun = None

  # Read image data from .npy file
  if kwargs['image_data_file'][-4:] == '.npy':
    image = np.load(kwargs['image_data_file'])

  # Read image data from .npz file
  elif kwargs['image_data_file'][-4:] == '.npz':
    data = np.load(kwargs['image_data_file'])
    image = data['image']
    if 'im_width' in data.keys():
      im_width = data['im_width'][0]
    if 'im_pos' in data.keys():
      if 'mass_msun' in data.keys():
        mass_msun = data['mass_msun'][0]

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

  # Calculate scale
  vmax = np.nanmax(image)
  scale_exponent = int('{0:24.16e}'.format(vmax).split('e')[1])
  scale = 10.0 ** scale_exponent
  image /= scale
  vmax /= scale
  if scale_exponent == 0:
    label = r'$I_\nu$ ($\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
  elif scale_exponent == 1:
    label = r'$I_\nu$ ($10\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'
  else:
    label = r'$I_\nu$ ($10^{' + repr(scale_exponent) \
        + r'}\ \mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}\ Hz^{-1}}$)'

  # Calculate grid
  extent = None
  if im_width is not None:
    extent = np.array((-im_width / 2.0, im_width / 2.0, -im_width / 2.0, im_width / 2.0))
    if mass_msun is not None and kwargs['distance'] is not None:
      r_g = gg_msun * mass_msun / c ** 2
      muas = np.pi / 180.0 / 60.0 / 60.0 / 1.0e6
      extent = np.arctan(extent * r_g / (kwargs['distance'] * pc)) / muas
      x_label = r'$x$ ($\mathrm{\mu as}$)'
      y_label = r'$y$ ($\mathrm{\mu as}$)'
    else:
      x_label = r'$x$ ($GM/c^2$)'
      y_label = r'$y$ ($GM/c^2$)'
  else:
    x_label = r'$x$-pixel'
    y_label = r'$y$-pixel'

  # Plot figure
  plt.imshow(image, cmap=cmap, vmin=0.0, vmax=vmax, aspect='equal', origin='lower', extent=extent,
      interpolation=interpolation)
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  plt.colorbar(label=label)
  plt.tight_layout()
  plt.savefig(kwargs['image_plot_file'], dpi=dpi)

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('image_data_file', help='file containing raw image data')
  parser.add_argument('image_plot_file', help='image file to create')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in parsecs')
  args = parser.parse_args()
  main(**vars(args))
