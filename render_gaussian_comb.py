#! /usr/bin/env python

"""
Script for creating rendering input parameters for Gaussian comb.
"""

# Python standard modules
import argparse
import warnings

# Numerical modules
import numpy as np
from scipy.optimize import bisect, least_squares
from scipy.stats import norm

# Plotting modules
import colorspacious as cs
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt

# Main function
def main(**kwargs):

  # Extract inputs
  quantity = kwargs['quantity']
  quantity_min = kwargs['quantity_min']
  quantity_max = kwargs['quantity_max']
  tau_scale_min = kwargs['tau_scale_min']
  tau_scale_max = kwargs['tau_scale_max']
  num_teeth = kwargs['num_teeth']
  rel_width = kwargs['rel_width']
  black_frac = kwargs['black_frac']
  num_tophats_color = kwargs['num_tophats_color']
  num_tophats_black = kwargs['num_tophats_black']
  sigma_max = kwargs['sigma_max']
  log = kwargs['log']
  log_scale = kwargs['log_scale']
  cmap = kwargs['cmap']
  offset = kwargs['offset']

  # Calculate range of quantity
  if log:
    quantity_min = np.log(quantity_min)
    quantity_max = np.log(quantity_max)

  # Calculate tooth centers
  centers, separation = np.linspace(quantity_min, quantity_max, num_teeth, retstep=True)

  # Calculate tooth width
  sigma = rel_width * separation

  # Calculate optical depths
  if log_scale:
    scales = np.exp(np.linspace(np.log(tau_scale_min), np.log(tau_scale_max), num_teeth))
  else:
    scales = np.linspace(tau_scale_min, tau_scale_max, num_teeth)

  # Calculate colors
  if cmap == 'viridis_alt':
    viridis_alt()
  colors = plt.get_cmap(cmap)(np.linspace(0.0, 1.0, num_teeth))

  # Calculate number of tophats
  if black_frac > 0.0:
    num_tophats_black = num_tophats_black
  else:
    num_tophats_black = 0
  num_tophats = num_tophats_color + 2 * num_tophats_black

  # Calculate tophat boundaries
  area = norm.cdf(sigma_max) - norm.cdf(-sigma_max)
  if black_frac > 0.0:
    sigmas_black = []
    sigmas_black.append(bisect(lambda x: norm.cdf(x) - 0.5 * (1.0 + area - black_frac * area), 0.0, sigma_max))
    for n in range(num_tophats_black - 1):
      sigmas_black.append(bisect(lambda x: norm.cdf(x) - 0.5 * (1.0 + area - black_frac * area + (n + 1.0) / num_tophats_black * black_frac * area), sigmas_black[0], sigma_max))
    sigmas_black.append(sigma_max)
  else:
    sigmas_black = [sigma_max]
  sigmas_color = []
  sigmas_color.append(-sigmas_black[0])
  for n in range(num_tophats_color - 1):
    sigmas_color.append(bisect(lambda x: norm.cdf(x) - 0.5 * (1.0 - area * (1.0 - black_frac) * (1.0 - 2.0 * (n + 1.0) / num_tophats_color)), sigmas_color[0], sigmas_black[0]))
  sigmas_color.append(sigmas_black[0])

  # Calculate tophat heights
  heights_color = []
  for n in range(num_tophats_color):
    heights_color.append((1.0 - black_frac) / num_tophats_color / (sigmas_color[n+1] - sigmas_color[n]))
  if black_frac > 0.0:
    heights_black = []
    for n in range(num_tophats_black):
      heights_black.append(0.5 * black_frac / num_tophats_color / (sigmas_black[n+1] - sigmas_black[n]))

  # Prepare output
  output_left = []
  output_right = []

  # Append initial output
  output_left.append('render_num_images')
  output_right.append('1')
  output_left.append('render_1_num_features')
  output_right.append(str(num_teeth * num_tophats))

  # Calculate black string
  black_string = '0.0,0.0,0.0'

  # Go through teeth
  for nn in range(num_teeth):

    # Calculate color string
    color_string = repr(colors[nn][0] * 255.0) + ',' + repr(colors[nn][1] * 255.0) + ',' + repr(colors[nn][2] * 255.0)

    # Go through left black tophats
    for n in range(num_tophats_black):

      # Calculate entry number
      n_entry = nn * num_tophats + n + offset + 1

      # Calculate range
      range_min = centers[nn] - sigmas_black[::-1][n] * sigma
      range_max = centers[nn] - sigmas_black[::-1][n+1] * sigma
      if log:
        range_min = np.exp(range_min)
        range_max = np.exp(range_max)

      # Calculate optical depth
      scale = scales[nn] / heights_black[n]

      # Append rendering entry to output
      output_left.append('render_1_' + str(n_entry) + '_type')
      output_right.append('fill')
      output_left.append('render_1_' + str(n_entry) + '_quantity')
      output_right.append(quantity)
      output_left.append('render_1_' + str(n_entry) + '_min')
      output_right.append(repr(range_min))
      output_left.append('render_1_' + str(n_entry) + '_max')
      output_right.append(repr(range_max))
      output_left.append('render_1_' + str(n_entry) + '_tau_scale')
      output_right.append(repr(scale))
      output_left.append('render_1_' + str(n_entry) + '_rgb')
      output_right.append(black_string)

    # Go through color tophats
    for n in range(num_tophats_color):

      # Calculate entry number
      n_entry = nn * num_tophats + num_tophats_black + n + offset + 1

      # Calculate range
      range_min = centers[nn] + sigmas_color[n] * sigma
      range_max = centers[nn] + sigmas_color[n+1] * sigma
      if log:
        range_min = np.exp(range_min)
        range_max = np.exp(range_max)

      # Calculate optical depth
      scale = scales[nn] / heights_color[n]

      # Append rendering entry to output
      output_left.append('render_1_' + str(n_entry) + '_type')
      output_right.append('fill')
      output_left.append('render_1_' + str(n_entry) + '_quantity')
      output_right.append(quantity)
      output_left.append('render_1_' + str(n_entry) + '_min')
      output_right.append(repr(range_min))
      output_left.append('render_1_' + str(n_entry) + '_max')
      output_right.append(repr(range_max))
      output_left.append('render_1_' + str(n_entry) + '_tau_scale')
      output_right.append(repr(scale))
      output_left.append('render_1_' + str(n_entry) + '_rgb')
      output_right.append(color_string)

    # Go through right black tophats
    for n in range(num_tophats_black):

      # Calculate entry number
      n_entry = nn * num_tophats + num_tophats_black + num_tophats_color + n + offset + 1

      # Calculate range
      range_min = centers[nn] + sigmas_black[n] * sigma
      range_max = centers[nn] + sigmas_black[n+1] * sigma
      if log:
        range_min = np.exp(range_min)
        range_max = np.exp(range_max)

      # Calculate optical depth
      scale = scales[nn] / heights_black[n]

      # Append rendering entry to output
      output_left.append('render_1_' + str(n_entry) + '_type')
      output_right.append('fill')
      output_left.append('render_1_' + str(n_entry) + '_quantity')
      output_right.append(quantity)
      output_left.append('render_1_' + str(n_entry) + '_min')
      output_right.append(repr(range_min))
      output_left.append('render_1_' + str(n_entry) + '_max')
      output_right.append(repr(range_max))
      output_left.append('render_1_' + str(n_entry) + '_tau_scale')
      output_right.append(repr(scale))
      output_left.append('render_1_' + str(n_entry) + '_rgb')
      output_right.append(black_string)

  # Print output
  print('')
  max_chars = np.max([len(output) for output in output_left]) + 1
  for left, right in zip(output_left, output_right):
    space = ' ' * (max_chars - len(left))
    print(left + space + '= ' + right)
  print('')

# Perceptually uniform monotonic colormap matching endpoints of viridis
def viridis_alt(name='viridis_alt', **kwargs):

  # Parameters
  cmap_base_name = 'viridis'
  start_val = 0.0
  end_val = 1.0
  winding_num = 0

  # Calculate endpoints
  start_rgb1 = cm.get_cmap(cmap_base_name)(start_val)[:-1]
  end_rgb1 = cm.get_cmap(cmap_base_name)(end_val)[:-1]
  start_jjmmh = cs.cspace_convert(start_rgb1, 'sRGB1', 'JMh')
  end_jjmmh = cs.cspace_convert(end_rgb1, 'sRGB1', 'JMh')

  # Create colormap
  return helix_uniform(start_jjmmh, end_jjmmh, winding_num, name=name, **kwargs)

# Perceptually uniform helix
def helix_uniform(start_jjmmh, end_jjmmh, winding_num, name='helix_uniform', num_points=1024):

  # Parameters
  c1 = 0.007
  c2 = 0.0228
  hr_max_dev = np.pi/4.0
  end_frac = 0.05

  # Prepare abscissas
  x = np.linspace(0.0, 1.0, num_points)

  # Calculate endpoints in perceptually uniform space
  jj_start = start_jjmmh[0]
  jjp_start = (1.0 + 100.0 * c1) * jj_start / (1.0 + c1 * jj_start)
  mm_start = start_jjmmh[1]
  mmp_start = np.log(1.0 + c2 * mm_start) / c2
  h_start = start_jjmmh[2]
  hr_start = h_start * np.pi / 180.0
  jj_end = end_jjmmh[0]
  jjp_end = (1.0 + 100.0 * c1) * jj_end / (1.0 + c1 * jj_end)
  mm_end = end_jjmmh[1]
  mmp_end = np.log(1.0 + c2 * mm_end) / c2
  h_end = end_jjmmh[2]
  hr_end = h_end * np.pi / 180.0

  # Calculate regular helix
  jjp = np.linspace(jjp_start, jjp_end, num_points)
  mmp = np.linspace(mmp_start, mmp_end, num_points)
  hr = np.linspace(hr_start, hr_end + 2.0 * np.pi * winding_num, num_points)

  # Make helix uniform in color difference
  def res(hr_int_vals):
    hr_vals = np.concatenate(([hr[0]], hr_int_vals, [hr[-1]]))
    ap_vals = mmp * np.cos(hr_vals)
    bp_vals = mmp * np.sin(hr_vals)
    deep_vals = np.hypot(np.diff(jjp), np.hypot(np.diff(ap_vals), np.diff(bp_vals)))
    return deep_vals - np.mean(deep_vals)
  def jac(hr_int_vals):
    hr_vals = np.concatenate(([hr[0]], hr_int_vals, [hr[-1]]))
    ap_vals = mmp * np.cos(hr_vals)
    bp_vals = mmp * np.sin(hr_vals)
    deep_vals = np.hypot(np.diff(jjp), np.hypot(np.diff(ap_vals), np.diff(bp_vals)))
    jacobian = np.zeros((num_points - 1, num_points - 2))
    jacobian[0,0] = mmp[0] * mmp[1] / deep_vals[0] * np.sin(hr_vals[1] - hr_vals[0])
    for k in range(1, num_points - 2):
      jacobian[k,k-1] = mmp[k] * mmp[k+1] / deep_vals[k] * np.sin(hr_vals[k] - hr_vals[k+1])
      jacobian[k,k] = -jacobian[k,k-1]
    jacobian[-1,-1] = mmp[-2] * mmp[-1] / deep_vals[-1] * np.sin(hr_vals[-1] - hr_vals[-2])
    return jacobian
  jac_sparsity = np.zeros((num_points - 1, num_points - 2))
  jac_sparsity[0,0] = 1.0
  for k in range(1, num_points - 2):
    jac_sparsity[k,k-1:k+1] = 1.0
  jac_sparsity[-1,-1] = 1.0
  bounds = (hr[1:-1] - hr_max_dev, hr[1:-1] + hr_max_dev)
  sol = least_squares(res, hr[1:-1], jac=jac, bounds=bounds, jac_sparsity=jac_sparsity)
  hr = np.concatenate(([hr[0]], sol.x, [hr[-1]]))

  # Calculate RGB values
  ap = mmp * np.cos(hr)
  bp = mmp * np.sin(hr)
  jjab = np.hstack((jjp[:,None], ap[:,None], bp[:,None]))
  with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'divide by zero encountered in divide', RuntimeWarning)
    rgb1 = cs.cspace_convert(jjab, cs.CAM02UCS, 'sRGB1')

  # Clip colors
  rgb1 = np.clip(rgb1, 0.0, 1.0)
  jjab = cs.cspace_convert(rgb1, 'sRGB1', cs.CAM02UCS)

  # Ensure clipped lightness is monotonic near ends
  ind_max = int(end_frac * (num_points - 1)) + 1
  if jjp[-1] - jjp[0] > 0.0 and np.any(np.diff(jjab[:ind_max,0]) < 0.0):
    ind = np.where(np.diff(jjab[:ind_max,0]) < 0.0)[0][-1] + 2
    fracs = (jjp[1:ind] - jjp[0]) / (jjp[ind] - jjp[0])
    rgb1[1:ind,:] = (1.0 - fracs[:,None]) * rgb1[:1,:] + fracs[:,None] * rgb1[ind:ind+1,:]
  elif jjp[-1] - jjp[0] < 0.0 and np.any(np.diff(jjab[:ind_max,0]) > 0.0):
    ind = np.where(np.diff(jjab[:ind_max,0]) > 0.0)[0][-1] + 2
    fracs = (jjp[1:ind] - jjp[0]) / (jjp[ind] - jjp[0])
    rgb1[1:ind,:] = (1.0 - fracs[:,None]) * rgb1[:1,:] + fracs[:,None] * rgb1[ind:ind+1,:]
  ind_min = int((1.0 - end_frac) * (num_points - 1))
  if jjp[-1] - jjp[0] > 0.0 and np.any(np.diff(jjab[ind_min:,0]) < 0.0):
    ind = np.where(np.diff(jjab[ind_min:,0]) < 0.0)[0][0] + ind_min - 1
    fracs = (jjp[ind+1:-1] - jjp[ind]) / (jjp[-1] - jjp[ind])
    rgb1[ind+1:-1,:] = (1.0 - fracs[:,None]) * rgb1[ind:ind+1,:] + fracs[:,None] * rgb1[-1:,:]
  elif jjp[-1] - jjp[0] < 0.0 and np.any(np.diff(jjab[ind_min:,0]) > 0.0):
    ind = np.where(np.diff(jjab[ind_min:,0]) > 0.0)[0][0] + ind_min - 1
    fracs = (jjp[ind+1:-1] - jjp[ind]) / (jjp[-1] - jjp[ind])
    rgb1[ind+1:-1,:] = (1.0 - fracs[:,None]) * rgb1[ind:ind+1,:] + fracs[:,None] * rgb1[-1:,:]

  # Create colormap
  x_rgb1 = [(x_val, rgb1_val) for x_val, rgb1_val in zip(x, rgb1)]
  cmap = colors.LinearSegmentedColormap.from_list(name, x_rgb1, N=num_points)
  cm.register_cmap(name=name, cmap=cmap)
  return cmap

# Parse inputs and execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('quantity', help='name of quantity to track')
  parser.add_argument('quantity_min', type=float, help='central value of lowest tooth')
  parser.add_argument('quantity_max', type=float, help='central value of highest tooth')
  parser.add_argument('tau_scale_min', type=float, help='optical depth scale of lowest tooth')
  parser.add_argument('tau_scale_max', type=float, help='optical depth scale of highest tooth')
  parser.add_argument('--num_teeth', type=int, default=5, help='number of teeth')
  parser.add_argument('--rel_width', type=float, default=0.1, help='sigma of teeth relative to tooth separation')
  parser.add_argument('--black_frac', type=float, default=0.5, help='fraction of Gaussian area set to black')
  parser.add_argument('--num_tophats_color', type=int, default=6, help='number of tophats to use in approximating central part of Gaussian')
  parser.add_argument('--num_tophats_black', type=int, default=3, help='number of tophats to use in approximating each wing of Gaussian')
  parser.add_argument('--sigma_max', type=float, default=3.0, help='maximum number of standard deviations to consider')
  parser.add_argument('--log', action='store_true', help='flag indicating comb should be constructed in log space')
  parser.add_argument('--log_scale', action='store_true', help='flag indicating optical depths should be constructed in log space')
  parser.add_argument('--cmap', default='viridis', help='name of colormap to use')
  parser.add_argument('--offset', type=int, default=0, help='offset to use for render features')
  args = parser.parse_args()
  main(**vars(args))
