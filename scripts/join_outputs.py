#! /usr/bin/env python

"""
Script for merging select datasets produced by Blacklight.
"""

# Python standard modules
import argparse

# Numerical modules
import numpy as np

# Main function
def main(**kwargs):

  # Verify inputs
  if kwargs['inputs'] is None:
    raise RuntimeError('Must supply at least one input file.')
  for filename_in in kwargs['inputs']:
    if filename_in[-4:] != '.npz':
      raise RuntimeError('Only npz format supported.')
  if kwargs['output'] is None:
    raise RuntimeError('Must supply name of output file.')
  if kwargs['output'][-4:] != '.npz':
    raise RuntimeError('Only npz format supported.')
  if kwargs['names'] is None:
    names_specified = False
    quantities = ('rho', 'n_e', 'p_gas', 'Theta_e', 'B', 'sigma', 'beta_inverse')
    names = ['positions', 'directions', 'I_nu', 'Q_nu', 'U_nu', 'V_nu', 'time', 'length', 'lambda',
        'emission', 'tau', 'rendering']
    names += ['lambda_ave_' + quantity for quantity in quantities]
    names += ['emission_ave_' + quantity for quantity in quantities]
    names += ['tau_int_' + quantity for quantity in quantities]
  else:
    names = kwargs['names']
    names_specified = True

  # Read non-adaptive data from first file
  data = {}
  with np.load(kwargs['inputs'][0], 'r') as f_in:
    data['mass_msun'] = f_in['mass_msun']
    data['width'] = f_in['width']
    data['frequency'] = f_in['frequency']
    names_present = []
    for name in names:
      try:
        data[name] = f_in[name]
      except KeyError as error:
        if names_specified:
          raise error
      else:
        names_present.append(name)
        if name == 'rendering':
          num_renderings = data[name].shape[0]

  # Extract metadata from all files
  num_freq = len(data['frequency'])
  multifreq = num_freq > 1
  adaptive_num_levels_all = {}
  adaptive_num_blocks_all = {}
  adaptive_block_locs_all = {}
  adaptive = False
  adaptive_block_size = None
  for filename_in in kwargs['inputs']:
    with np.load(filename_in, 'r') as f_in:
      adaptive_num_levels_all[filename_in] = f_in['adaptive_num_levels']
      try:
        adaptive_num_blocks_all[filename_in] = f_in['adaptive_num_blocks']
      except KeyError:
        adaptive_num_blocks_all[filename_in] = None
      else:
        adaptive = True
      adaptive_block_locs_all[filename_in] = {}
      for n in range(1, adaptive_num_levels_all[filename_in] + 1):
        adaptive_block_locs_all[filename_in][n] = f_in['adaptive_block_locs_{0}'.format(n)]
        if adaptive_block_size is None:
          if names_present[0] in ('positions', 'directions', 'time', 'length'):
            adaptive_block_size = f_in['adaptive_{0}_{1}'.format(names_present[0],n)].shape[1]
          elif names_present[0] == 'rendering':
            adaptive_block_size = f_in['adaptive_{0}_{1}'.format(names_present[0],n)].shape[3]
          else:
            if multifreq:
              adaptive_block_size = f_in['adaptive_{0}_{1}'.format(names_present[0],n)].shape[2]
            else:
              adaptive_block_size = f_in['adaptive_{0}_{1}'.format(names_present[0],n)].shape[1]

  # Calculate maximum refinement across all files
  adaptive_num_levels = 0
  for filename_in in kwargs['inputs']:
    adaptive_num_levels = max(adaptive_num_levels, adaptive_num_levels_all[filename_in])

  # Prepare description of joined data
  adaptive_num_blocks = np.zeros(adaptive_num_levels + 1, dtype=int)
  adaptive_block_locs = {}
  for n in range(1, adaptive_num_levels + 1):
    adaptive_block_locs[n] = []

  # Determine which blocks need extracting
  block_flags_all = {}
  for filename_in in kwargs['inputs']:
    block_flags_all[filename_in] = {}
    for n in range(1, adaptive_num_levels_all[filename_in] + 1):
      block_flags_all[filename_in][n] = []
      for block in range(adaptive_num_blocks_all[filename_in][n]):
        yx_loc = list(adaptive_block_locs_all[filename_in][n][block,:])
        if yx_loc not in adaptive_block_locs[n]:
          block_flags_all[filename_in][n].append(block)
          adaptive_num_blocks[n] += 1
          adaptive_block_locs[n].append(yx_loc)

  # Store adaptive metadata
  adaptive_num_levels_array = np.empty(1, dtype=int)
  adaptive_num_levels_array[0] = adaptive_num_levels
  data['adaptive_num_levels'] = adaptive_num_levels_array
  if adaptive:
    data['adaptive_num_blocks'] = adaptive_num_blocks
    for n in range(1, adaptive_num_levels + 1):
      data['adaptive_block_locs_{0}'.format(n)] = adaptive_block_locs[n]

  # Prepare joined adaptive data
  for n in range(1, adaptive_num_levels + 1):
    for name in names_present:
      if name in ('positions', 'directions'):
        shape = (adaptive_num_blocks[n], adaptive_block_size, adaptive_block_size, 4)
      elif name in ('time', 'length'):
        shape = (adaptive_num_blocks[n], adaptive_block_size, adaptive_block_size)
      elif name == 'rendering':
        shape = \
            (num_renderings, 3, adaptive_num_blocks[n], adaptive_block_size, adaptive_block_size)
      else:
        if multifreq:
          shape = (num_freq, adaptive_num_blocks[n], adaptive_block_size, adaptive_block_size)
        else:
          shape = (adaptive_num_blocks[n], adaptive_block_size, adaptive_block_size)
      data['adaptive_{0}_{1}'.format(name,n)] = np.empty(shape)

  # Read non-duplicate adaptive data from all files
  block_counts = np.zeros(adaptive_num_levels + 1, dtype=int)
  for filename_in in kwargs['inputs']:
    with np.load(filename_in, 'r') as f_in:
      for n in range(1, adaptive_num_levels_all[filename_in] + 1):
        for block in block_flags_all[filename_in][n]:
          for name in names_present:
            name_local = 'adaptive_{0}_{1}'.format(name, n)
            if name in ('positions', 'directions'):
              data[name_local][block_counts[n],:,:,:] = f_in[name_local][block,:,:,:]
            elif name in ('time', 'length'):
              data[name_local][block_counts[n],:,:] = f_in[name_local][block,:,:]
            elif name == 'rendering':
              data[name_local][:,:,block_counts[n],:,:] = f_in[name_local][:,:,block,:,:]
            else:
              if multifreq:
                data[name_local][:,block_counts[n],:,:] = f_in[name_local][:,block,:,:]
              else:
                data[name_local][block_counts[n],:,:] = f_in[name_local][block,:,:]
          block_counts[n] += 1

  # Save data
  np.savez(kwargs['output'], **data)

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--inputs', nargs='+', help='names of input files to read')
  parser.add_argument('-o', '--output', help='name of output file to write')
  parser.add_argument('-n', '--names', nargs='+', help='names of quantities to extract')
  args = parser.parse_args()
  main(**vars(args))
