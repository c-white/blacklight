#! /usr/bin/env python

"""
Script for plotting outputs produced by Blacklight as true-color images.
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

  # Verify input
  if kwargs['filename_data'][-4:] != '.npz':
    raise RuntimeError('Only npz format supports multiple frequencies.')
  if kwargs['central_frequency'] is None:
    raise RuntimeError('Must specify central_frequency.')

  # Prepare metadata
  distance_pc = kwargs['distance']
  max_level = kwargs['max_level']

  # Read data from .npz file
  with np.load(kwargs['filename_data']) as f:

    # Read metadata
    width_rg = f['width'][0]
    mass_msun = f['mass_msun'][0]
    frequencies = f['frequency'][:]

    # Read root intensity
    try:
      intensity = f['I_nu'][:]
    except KeyError:
      raise RuntimeError('No intensity data in file.')

    # Read adaptive intensity
    if max_level is None:
      max_level = f['adaptive_num_levels'][0]
    else:
      max_level = min(max_level, f['adaptive_num_levels'][0])
    if max_level > 0:
      num_blocks = {}
      block_locs = {}
      intensity_adaptive = {}
      num_blocks[0] = f['adaptive_num_blocks'][0]
      for level in range(1, max_level + 1):
        num_blocks[level] = f['adaptive_num_blocks'][level]
        block_locs[level] = f['adaptive_block_locs_{0}'.format(level)][:]
        intensity_adaptive[level] = f['adaptive_I_nu_{0}'.format(level)][:]

  # Calculate shifted wavelengths
  wavelengths_nm = kwargs['central_frequency'] * kwargs['central_wavelength'] / frequencies

  # Ensure data is ordered by increasing wavelength
  if np.all(np.diff(wavelengths_nm) > 0.0):
    pass
  elif np.all(np.diff(wavelengths_nm) < 0.0):
    wavelengths_nm = wavelengths_nm[::-1]
    intensity = intensity[::-1,...]
    for level in range(1, max_level + 1):
      intensity_adaptive[level] = intensity_adaptive[level][::-1,...]
  else:
    raise RuntimeError('Non-monotonic wavelengths not supported.')

  # Transform I_nu to I_lambda
  intensity /= wavelengths_nm[:,None,None] ** 2
  for level in range(1, max_level + 1):
    intensity_adaptive[level] /= wavelengths_nm[:,None,None,None] ** 2

  # Calculate RGB colors
  xyz = intensity_to_xyz(wavelengths_nm, intensity, True)
  rgb = xyz_to_rgb(xyz)
  image = np.concatenate((rgb[0,:,:,None], rgb[1,:,:,None], rgb[2,:,:,None]), axis=2)
  if max_level > 0:
    image_adaptive = {}
    for level in range(1, max_level + 1):
      xyz = intensity_to_xyz(wavelengths_nm, intensity_adaptive[level], False)
      rgb = xyz_to_rgb(xyz)
      image_adaptive[level] = \
          np.concatenate((rgb[0,:,:,:,None], rgb[1,:,:,:,None], rgb[2,:,:,:,None]), axis=3)

  # Calculate root grid in pixels
  if kwargs['axes'] == 'pixel':
    extent = np.array((-0.5, image.shape[0] - 0.5, -0.5, image.shape[1] - 0.5))
    x_label = r'$x$-pixel'
    y_label = r'$y$-pixel'

  # Calculate root grid in gravitational radii
  elif kwargs['axes'] == 'rg':
    if width_rg is None:
      raise RuntimeError('Must supply width.')
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
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    if mass_msun is None:
      raise RuntimeError('Must supply mass.')
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
    if width_rg is None:
      raise RuntimeError('Must supply width.')
    if mass_msun is None:
      raise RuntimeError('Must supply mass.')
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

# Function for converting intensities to CIE XYZ colors
# Notes:
#   Uses CIE 170-2 (2-degree) rather than CIE 1931 definitions:
#     L, M, and S cone response functions are from 2000 Vision Research 40 1711.
#     Y matching function is luminous efficiency from 2005 Journal of Vision 5 948 and 2011 Color
#         Research and Application 36 42.
#     Z matching function is scaled version of S cone response function from 1999 Vision Research 39
#         2901.
#     X matching function is chosen to minimize differences with CIE 1931.
def intensity_to_xyz(wavelengths, intensities, warn):

  # Define wavelengths
  match_lambdas = np.linspace(390.0, 830.0, 441)

  # Verify inputs
  if warn:
    if wavelengths[0] > match_lambdas[0]:
      deficit = wavelengths[0] - match_lambdas[0]
      if deficit > min(0.5 * (wavelengths[1] - wavelengths[0]), 25.0):
        warning_1 = 'Warning: shortest wavelength ({0} nm)'.format(wavelengths[0])
        warning_2 = ' significantly greater than expected ({0} nm)'.format(match_lambdas[0])
        print(warning_1 + warning_2)
    if wavelengths[-1] < match_lambdas[-1]:
      deficit = match_lambdas[-1] - wavelengths[-1]
      if deficit > min(0.5 * (wavelengths[-1] - wavelengths[-2]), 25.0):
        warning_1 = 'Warning: longest wavelength ({0} nm)'.format(wavelengths[-1])
        warning_2 = ' significantly less than expected ({0} nm)'.format(match_lambdas[-1])
        print(warning_1 + warning_2)

  # Define cone response functions
  cone_l = np.array((4.15003360e-4, 5.02650370e-4, 6.07366773e-4, 7.31850005e-4, 8.79011508e-4,
      1.05192281e-3, 1.25373403e-3, 1.48756033e-3, 1.75633272e-3, 2.06261148e-3, 2.40836274e-3,
      2.79522061e-3, 3.22639729e-3, 3.70617133e-3, 4.23972185e-3, 4.83338951e-3, 5.49335058e-3,
      6.21932550e-3, 7.00630870e-3, 7.84503388e-3, 8.72126562e-3, 9.61879153e-3, 1.05323536e-2,
      1.14619761e-2, 1.24104616e-2, 1.33836930e-2, 1.43870208e-2, 1.54116155e-2, 1.64424078e-2,
      1.74614455e-2, 1.84480328e-2, 1.93852304e-2, 2.02811387e-2, 2.11544610e-2, 2.20286383e-2,
      2.29317499e-2, 2.38896150e-2, 2.49026216e-2, 2.59631430e-2, 2.70618851e-2, 2.81877052e-2,
      2.93303393e-2, 3.04898213e-2, 3.16694350e-2, 3.28730969e-2, 3.41053893e-2, 3.53670491e-2,
      3.66400331e-2, 3.78988484e-2, 3.91148083e-2, 4.02562524e-2, 4.12989224e-2, 4.22582196e-2,
      4.31627409e-2, 4.40443720e-2, 4.49379920e-2, 4.58728905e-2, 4.68459041e-2, 4.78446315e-2,
      4.88555461e-2, 4.98639402e-2, 5.08618457e-2, 5.18730928e-2, 5.29317333e-2, 5.40746201e-2,
      5.53418498e-2, 5.67734433e-2, 5.83973435e-2, 6.02408731e-2, 6.23352567e-2, 6.47164360e-2,
      6.74130612e-2, 7.04040119e-2, 7.36488902e-2, 7.70978002e-2, 8.06893625e-2, 8.43613293e-2,
      8.80904420e-2, 9.18629938e-2, 9.56636658e-2, 9.94755270e-2, 1.03285638e-1, 1.07103485e-1,
      1.10947170e-1, 1.14838313e-1, 1.18802152e-1, 1.22862603e-1, 1.27025792e-1, 1.31293061e-1,
      1.35665739e-1, 1.40145131e-1, 1.44730863e-1, 1.49415456e-1, 1.54188534e-1, 1.59038311e-1,
      1.63951537e-1, 1.68932334e-1, 1.74063331e-1, 1.79457342e-1, 1.85240873e-1, 1.91556224e-1,
      1.98532355e-1, 2.06182424e-1, 2.14484817e-1, 2.23411375e-1, 2.32925516e-1, 2.42992175e-1,
      2.53615531e-1, 2.64809627e-1, 2.76587210e-1, 2.88959412e-1, 3.01933704e-1, 3.15507769e-1,
      3.29673082e-1, 3.44415858e-1, 3.59716421e-1, 3.75550447e-1, 3.91894602e-1, 4.08721523e-1,
      4.25997662e-1, 4.43682863e-1, 4.61731088e-1, 4.80093753e-1, 4.98716765e-1, 5.17539319e-1,
      5.36493787e-1, 5.55492678e-1, 5.74386124e-1, 5.92994956e-1, 6.11123546e-1, 6.28561087e-1,
      6.45138591e-1, 6.60907396e-1, 6.75993265e-1, 6.90542475e-1, 7.04720332e-1, 7.18662310e-1,
      7.32322897e-1, 7.45605717e-1, 7.58409759e-1, 7.70629906e-1, 7.82207895e-1, 7.93289637e-1,
      8.04083976e-1, 8.14812845e-1, 8.25710771e-1, 8.36942718e-1, 8.48348797e-1, 8.59675605e-1,
      8.70656228e-1, 8.81010908e-1, 8.90495937e-1, 8.99051792e-1, 9.06668909e-1, 9.13341335e-1,
      9.19066658e-1, 9.23897675e-1, 9.28099312e-1, 9.31995493e-1, 9.35916341e-1, 9.40197682e-1,
      9.45075812e-1, 9.50366568e-1, 9.55775099e-1, 9.60999943e-1, 9.65733133e-1, 9.69744189e-1,
      9.73133076e-1, 9.76086087e-1, 9.78792641e-1, 9.81444621e-1, 9.84186485e-1, 9.86965393e-1,
      9.89678050e-1, 9.92220033e-1, 9.94485919e-1, 9.96385799e-1, 9.97894643e-1, 9.99004045e-1,
      9.99706046e-1, 9.99993160e-1, 9.99836906e-1, 9.99123467e-1, 9.97718715e-1, 9.95490549e-1,
      9.92309909e-1, 9.88145905e-1, 9.83345096e-1, 9.78342318e-1, 9.73564193e-1, 9.69429172e-1,
      9.66211047e-1, 9.63629933e-1, 9.61272895e-1, 9.58731517e-1, 9.55601819e-1, 9.51568637e-1,
      9.46653575e-1, 9.40962330e-1, 9.34600447e-1, 9.27672559e-1, 9.20264315e-1, 9.12390740e-1,
      9.04050350e-1, 8.95242882e-1, 8.85969334e-1, 8.76242466e-1, 8.66117284e-1, 8.55657805e-1,
      8.44926107e-1, 8.33982132e-1, 8.22858885e-1, 8.11491117e-1, 7.99794023e-1, 7.87688857e-1,
      7.75103324e-1, 7.61995930e-1, 7.48424747e-1, 7.34469825e-1, 7.20208108e-1, 7.05713103e-1,
      6.91043503e-1, 6.76211520e-1, 6.61219505e-1, 6.46071535e-1, 6.30773388e-1, 6.15349220e-1,
      5.99888258e-1, 5.84488739e-1, 5.69240353e-1, 5.54224372e-1, 5.39469305e-1, 5.24827070e-1,
      5.10123556e-1, 4.95205531e-1, 4.79940995e-1, 4.64269520e-1, 4.48338311e-1, 4.32329322e-1,
      4.16405540e-1, 4.00710831e-1, 3.85355429e-1, 3.70376724e-1, 3.55792605e-1, 3.41617724e-1,
      3.27863707e-1, 3.14540835e-1, 3.01662231e-1, 2.89238816e-1, 2.77278035e-1, 2.65784185e-1,
      2.54739669e-1, 2.44054451e-1, 2.33633978e-1, 2.23398785e-1, 2.13283643e-1, 2.03256632e-1,
      1.93370014e-1, 1.83687512e-1, 1.74262972e-1, 1.65140826e-1, 1.56353734e-1, 1.47915648e-1,
      1.39833545e-1, 1.32111108e-1, 1.24749078e-1, 1.17743591e-1, 1.11081356e-1, 1.04747491e-1,
      9.87277483e-2, 9.30084863e-2, 8.75769371e-2, 8.24219445e-2, 7.75328030e-2, 7.28989439e-2,
      6.85099590e-2, 6.43554002e-2, 6.04242407e-2, 5.67057277e-2, 5.31896065e-2, 4.98660938e-2,
      4.67258325e-2, 4.37598140e-2, 4.09594391e-2, 3.83165189e-2, 3.58232509e-2, 3.34724721e-2,
      3.12583135e-2, 2.91751380e-2, 2.72172710e-2, 2.53790394e-2, 2.36540832e-2, 2.20336126e-2,
      2.05092130e-2, 1.90735277e-2, 1.77201275e-2, 1.64450960e-2, 1.52510378e-2, 1.41403568e-2,
      1.31136596e-2, 1.21701123e-2, 1.13063241e-2, 1.05131103e-2, 9.78126875e-3, 9.10300059e-3,
      8.47170037e-3, 7.88226553e-3, 7.33213030e-3, 6.81928276e-3, 6.34172502e-3, 5.89748808e-3,
      5.48443753e-3, 5.09977946e-3, 4.74086432e-3, 4.40537654e-3, 4.09129197e-3, 3.79702858e-3,
      3.52186853e-3, 3.26519932e-3, 3.02632064e-3, 2.80446659e-3, 2.59881597e-3, 2.40848695e-3,
      2.23258808e-3, 2.07023891e-3, 1.92057783e-3, 1.78268721e-3, 1.65540067e-3, 1.53761987e-3,
      1.42838689e-3, 1.32686544e-3, 1.23238346e-3, 1.14456181e-3, 1.06307564e-3, 9.87592383e-4,
      9.17776773e-4, 8.53264265e-4, 7.93589168e-4, 7.38305966e-4, 6.87018276e-4, 6.39372955e-4,
      5.95056902e-4, 5.53797304e-4, 5.15350300e-4, 4.79495834e-4, 4.46035067e-4, 4.14809228e-4,
      3.85749126e-4, 3.58791917e-4, 3.33860879e-4, 3.10869146e-4, 2.89698852e-4, 2.70144543e-4,
      2.52007194e-4, 2.35117002e-4, 2.19329240e-4, 2.04534887e-4, 1.90692352e-4, 1.77770930e-4,
      1.65735819e-4, 1.54549195e-4, 1.44166855e-4, 1.34528395e-4, 1.25574343e-4, 1.17250536e-4,
      1.09507604e-4, 1.02299795e-4, 9.55827969e-5, 8.93160676e-5, 8.34631462e-5, 7.79912028e-5,
      7.28729665e-5, 6.80921255e-5, 6.36341781e-5, 5.94840578e-5, 5.56263839e-5, 5.20429804e-5,
      4.87063495e-5, 4.55900001e-5, 4.26709878e-5, 3.99294560e-5, 3.73509304e-5, 3.49326465e-5,
      3.26730386e-5, 3.05690170e-5, 2.86162863e-5, 2.68076897e-5, 2.51286065e-5, 2.35644569e-5,
      2.21026115e-5, 2.07321327e-5, 1.94444103e-5, 1.82351096e-5, 1.71007926e-5, 1.60379855e-5,
      1.50432118e-5, 1.41127462e-5, 1.32419536e-5, 1.24263363e-5, 1.16618002e-5, 1.09446148e-5,
      1.02715618e-5, 9.64036490e-6, 9.04896926e-6, 8.49534579e-6, 7.97750298e-6, 7.49341023e-6,
      7.04078840e-6, 6.61743908e-6, 6.22132818e-6, 5.85057291e-6, 5.50333412e-6, 5.17756349e-6,
      4.87135044e-6, 4.58300170e-6, 4.31101652e-6, 4.05421999e-6, 3.81213497e-6, 3.58438418e-6,
      3.37052966e-6, 3.17008616e-6, 2.98247763e-6, 2.80690747e-6, 2.64257450e-6, 2.48873179e-6,
      2.34468299e-6, 2.20974413e-6, 2.08315199e-6, 1.96419005e-6, 1.85221800e-6, 1.74666331e-6,
      1.64704844e-6, 1.55307208e-6, 1.46447518e-6, 1.38100319e-6, 1.30240694e-6, 1.22844348e-6,
      1.15887663e-6, 1.09347757e-6, 1.03202519e-6, 9.74306403e-7))
  cone_m = np.array((3.68349248e-4, 4.48014903e-4, 5.43965499e-4, 6.58983162e-4, 7.96121210e-4,
      9.58658287e-4, 1.15002180e-3, 1.37367427e-3, 1.63295701e-3, 1.93088710e-3, 2.26990665e-3,
      2.65210441e-3, 3.08110117e-3, 3.56155706e-3, 4.09900407e-3, 4.70010083e-3, 5.37186059e-3,
      6.11757407e-3, 6.93795306e-3, 7.83144260e-3, 8.79368990e-3, 9.81865357e-3, 1.09044480e-2,
      1.20508507e-2, 1.32582230e-2, 1.45277327e-2, 1.58600151e-2, 1.72495960e-2, 1.86877573e-2,
      2.01637692e-2, 2.16648593e-2, 2.31800771e-2, 2.47139091e-2, 2.62784761e-2, 2.78904592e-2,
      2.95714420e-2, 3.13437635e-2, 3.32150920e-2, 3.51888593e-2, 3.72683413e-2, 3.94566147e-2,
      4.17545976e-2, 4.41543810e-2, 4.66432409e-2, 4.92050252e-2, 5.18198823e-2, 5.44644516e-2,
      5.71130840e-2, 5.97369425e-2, 6.23037817e-2, 6.47782187e-2, 6.71329339e-2, 6.93844290e-2,
      7.15658156e-2, 7.37163628e-2, 7.58812175e-2, 7.80979917e-2, 8.03537509e-2, 8.26198303e-2,
      8.48643801e-2, 8.70523965e-2, 8.91639889e-2, 9.12522896e-2, 9.33950398e-2, 9.56775066e-2,
      9.81933595e-2, 1.01031383e-1, 1.04229303e-1, 1.07814112e-1, 1.11816660e-1, 1.16272079e-1,
      1.21204305e-1, 1.26573030e-1, 1.32310541e-1, 1.38333533e-1, 1.44540604e-1, 1.50827804e-1,
      1.57145687e-1, 1.63457016e-1, 1.69720863e-1, 1.75892769e-1, 1.81937722e-1, 1.87871999e-1,
      1.93730862e-1, 1.99556529e-1, 2.05398073e-1, 2.11302131e-1, 2.17283004e-1, 2.23347021e-1,
      2.29501355e-1, 2.35754060e-1, 2.42108491e-1, 2.48545289e-1, 2.55037126e-1, 2.61553968e-1,
      2.68063038e-1, 2.74561420e-1, 2.81180052e-1, 2.88097133e-1, 2.95508395e-1, 3.03629653e-1,
      3.12649158e-1, 3.22564974e-1, 3.33319240e-1, 3.44844385e-1, 3.57060581e-1, 3.69895499e-1,
      3.83355071e-1, 3.97467012e-1, 4.12260159e-1, 4.27764487e-1, 4.44000827e-1, 4.60946806e-1,
      4.78561797e-1, 4.96795479e-1, 5.15586863e-1, 5.34867541e-1, 5.54574593e-1, 5.74639585e-1,
      5.94984257e-1, 6.15520174e-1, 6.36162030e-1, 6.56871646e-1, 6.77623604e-1, 6.98392699e-1,
      7.19154047e-1, 7.39844183e-1, 7.60235488e-1, 7.80038616e-1, 7.98940996e-1, 8.16609632e-1,
      8.32782761e-1, 8.47550216e-1, 8.61113469e-1, 8.73697924e-1, 8.85549719e-1, 8.96872538e-1,
      9.07638015e-1, 9.17756577e-1, 9.27137002e-1, 9.35687022e-1, 9.43361575e-1, 9.50308999e-1,
      9.56732838e-1, 9.62843704e-1, 9.68858338e-1, 9.74919067e-1, 9.80849502e-1, 9.86387994e-1,
      9.91267411e-1, 9.95216550e-1, 9.98007488e-1, 9.99594634e-1, 9.99981748e-1, 9.99176952e-1,
      9.97192616e-1, 9.94101614e-1, 9.90204438e-1, 9.85855446e-1, 9.81403969e-1, 9.77193335e-1,
      9.73440778e-1, 9.69880983e-1, 9.66133321e-1, 9.61823142e-1, 9.56583237e-1, 9.50166600e-1,
      9.42773123e-1, 9.34708698e-1, 9.26271394e-1, 9.17749837e-1, 9.09338868e-1, 9.00894806e-1,
      8.92197240e-1, 8.83034571e-1, 8.73204873e-1, 8.62565296e-1, 8.51172855e-1, 8.39131971e-1,
      8.26544532e-1, 8.13509255e-1, 8.00082372e-1, 7.86166191e-1, 7.71635006e-1, 7.56376136e-1,
      7.40291037e-1, 7.23368777e-1, 7.05889868e-1, 6.88184300e-1, 6.70554222e-1, 6.53273547e-1,
      6.36524204e-1, 6.20217289e-1, 6.04210507e-1, 5.88375391e-1, 5.72596792e-1, 5.56783148e-1,
      5.40895608e-1, 5.24912224e-1, 5.08816788e-1, 4.92598798e-1, 4.76267813e-1, 4.59892636e-1,
      4.43550966e-1, 4.27313848e-1, 4.11245617e-1, 3.95396356e-1, 3.79781665e-1, 3.64409629e-1,
      3.49288891e-1, 3.34428550e-1, 3.19842911e-1, 3.05564226e-1, 2.91625076e-1, 2.78053259e-1,
      2.64872009e-1, 2.52099488e-1, 2.39747055e-1, 2.27822335e-1, 2.16330273e-1, 2.05273346e-1,
      1.94649648e-1, 1.84447615e-1, 1.74654040e-1, 1.65256374e-1, 1.56242654e-1, 1.47602388e-1,
      1.39328896e-1, 1.31415765e-1, 1.23855877e-1, 1.16641491e-1, 1.09765763e-1, 1.03226432e-1,
      9.70205151e-2, 9.11430058e-2, 8.55871752e-2, 8.03427998e-2, 7.53912660e-2, 7.07135564e-2,
      6.62923873e-2, 6.21120321e-2, 5.81594780e-2, 5.44275517e-2, 5.09097902e-2, 4.75990951e-2,
      4.44878727e-2, 4.15659453e-2, 3.88151890e-2, 3.62180917e-2, 3.37598910e-2, 3.14282295e-2,
      2.92175652e-2, 2.71403387e-2, 2.52084169e-2, 2.34286294e-2, 2.18037465e-2, 2.03284511e-2,
      1.89779190e-2, 1.77272061e-2, 1.65560098e-2, 1.54479552e-2, 1.43923782e-2, 1.33895734e-2,
      1.24414195e-2, 1.15488447e-2, 1.07119605e-2, 9.92994435e-3, 9.20056706e-3, 8.52126391e-3,
      7.88945111e-3, 7.30255150e-3, 6.75820484e-3, 6.25473908e-3, 5.79045265e-3, 5.36346743e-3,
      4.97179190e-3, 4.61300300e-3, 4.28338721e-3, 3.97940401e-3, 3.69802741e-3, 3.43667061e-3,
      3.19326237e-3, 2.96653255e-3, 2.75542939e-3, 2.55895674e-3, 2.37617227e-3, 2.20617719e-3,
      2.04809376e-3, 1.90109614e-3, 1.76441584e-3, 1.63733769e-3, 1.51916415e-3, 1.40913442e-3,
      1.30654517e-3, 1.21077578e-3, 1.12127805e-3, 1.03765879e-3, 9.59891569e-4, 8.87951392e-4,
      8.21729099e-4, 7.61050801e-4, 7.05623212e-4, 6.54873950e-4, 6.08240457e-4, 5.65239773e-4,
      5.25456816e-4, 4.88549556e-4, 4.54277238e-4, 4.22435937e-4, 3.92839391e-4, 3.65317316e-4,
      3.39709748e-4, 3.15855564e-4, 2.93607253e-4, 2.72833505e-4, 2.53417295e-4, 2.35264961e-4,
      2.18330806e-4, 2.02574019e-4, 1.87947963e-4, 1.74401753e-4, 1.61877983e-4, 1.50305030e-4,
      1.39612176e-4, 1.29733576e-4, 1.20607982e-4, 1.12175830e-4, 1.04372509e-4, 9.71382564e-5,
      9.04202287e-5, 8.41716373e-5, 7.83538244e-5, 7.29425072e-5, 6.79164009e-5, 6.32542269e-5,
      5.89349310e-5, 5.49360120e-5, 5.12291125e-5, 4.77871885e-5, 4.45862207e-5, 4.16048509e-5,
      3.88245017e-5, 3.62301807e-5, 3.38085674e-5, 3.15474068e-5, 2.94354152e-5, 2.74633548e-5,
      2.56268159e-5, 2.39217211e-5, 2.23432310e-5, 2.08859570e-5, 1.95424279e-5, 1.82989666e-5,
      1.71422768e-5, 1.60610643e-5, 1.50457558e-5, 1.40893453e-5, 1.31898859e-5, 1.23462272e-5,
      1.15568653e-5, 1.08200177e-5, 1.01334139e-5, 9.49366644e-6, 8.89736081e-6, 8.34134639e-6,
      7.82271489e-6, 7.33865490e-6, 6.88611850e-6, 6.46227717e-6, 6.06462194e-6, 5.69092610e-6,
      5.33942055e-6, 5.00928587e-6, 4.69984414e-6, 4.41034426e-6, 4.13998115e-6, 3.88772428e-6,
      3.65185801e-6, 3.43069961e-6, 3.22278205e-6, 3.02682745e-6, 2.84191652e-6, 2.66795139e-6,
      2.50491424e-6, 2.35267302e-6, 2.21100363e-6, 2.07945815e-6, 1.95699850e-6, 1.84257587e-6,
      1.73528111e-6, 1.63432668e-6, 1.53909945e-6, 1.44932010e-6, 1.36477800e-6, 1.28525715e-6,
      1.21053851e-6, 1.14038010e-6, 1.07446298e-6, 1.01247449e-6, 9.54129935e-7, 8.99169870e-7,
      8.47378525e-7, 7.98634599e-7, 7.52831739e-7, 7.09857251e-7, 6.69594072e-7, 6.31906875e-7,
      5.96603670e-7, 5.63494696e-7, 5.32407946e-7, 5.03187359e-7, 4.75686244e-7, 4.49753130e-7,
      4.25248753e-7, 4.02049638e-7, 3.80046377e-7, 3.59156013e-7, 3.39356264e-7, 3.20632948e-7,
      3.02965616e-7, 2.86328659e-7, 2.70686503e-7, 2.55979580e-7, 2.42146542e-7, 2.29130415e-7,
      2.16878259e-7, 2.05337763e-7, 1.94448963e-7, 1.84155347e-7, 1.74406564e-7, 1.65157769e-7,
      1.56372504e-7, 1.48031083e-7, 1.40117420e-7, 1.32615288e-7, 1.25508416e-7, 1.18780584e-7,
      1.12415697e-7, 1.06397865e-7, 1.00711454e-7, 9.53411481e-8))
  cone_s = np.array((9.54728799e-3, 1.14794259e-2, 1.37985579e-2, 1.65745750e-2, 1.98869445e-2,
      2.38249757e-2, 2.84877010e-2, 3.39832301e-2, 4.04274189e-2, 4.79416817e-2, 5.66497751e-2,
      6.66756827e-2, 7.81479086e-2, 9.11924744e-2, 1.05926235e-1, 1.22450936e-1, 1.40843885e-1,
      1.61139855e-1, 1.83324677e-1, 2.07326960e-1, 2.33008318e-1, 2.60183226e-1, 2.88722693e-1,
      3.18511681e-1, 3.49431207e-1, 3.81363350e-1, 4.14140896e-1, 4.47350309e-1, 4.80439025e-1,
      5.12767135e-1, 5.43618105e-1, 5.72399437e-1, 5.99284259e-1, 6.24785633e-1, 6.49575739e-1,
      6.74474074e-1, 7.00186409e-1, 7.26460271e-1, 7.52725805e-1, 7.78333016e-1, 8.02554955e-1,
      8.24817858e-1, 8.45422394e-1, 8.64961201e-1, 8.84100149e-1, 9.03572547e-1, 9.23844418e-1,
      9.44055273e-1, 9.62919656e-1, 9.79057060e-1, 9.91020299e-1, 9.97764785e-1, 9.99981961e-1,
      9.98860947e-1, 9.95628403e-1, 9.91515250e-1, 9.87377323e-1, 9.82619286e-1, 9.76301170e-1,
      9.67512827e-1, 9.55393000e-1, 9.39498618e-1, 9.20807415e-1, 9.00592250e-1, 8.80040181e-1,
      8.60240453e-1, 8.42031293e-1, 8.25571506e-1, 8.10860447e-1, 7.97901650e-1, 7.86703893e-1,
      7.77118717e-1, 7.68364561e-1, 7.59537988e-1, 7.49777186e-1, 7.38267621e-1, 7.24389338e-1,
      7.08112614e-1, 6.89571829e-1, 6.68926552e-1, 6.46359219e-1, 6.22112167e-1, 5.96590777e-1,
      5.70215640e-1, 5.43373087e-1, 5.16410649e-1, 4.89647329e-1, 4.63405762e-1, 4.37965017e-1,
      4.13549483e-1, 3.90333407e-1, 3.68393733e-1, 3.47581782e-1, 3.27723819e-1, 3.08676241e-1,
      2.90322411e-1, 2.72626410e-1, 2.55771784e-1, 2.39945108e-1, 2.25280158e-1, 2.11866834e-1,
      1.99718535e-1, 1.88674025e-1, 1.78556779e-1, 1.69216691e-1, 1.60525792e-1, 1.52368238e-1,
      1.44620438e-1, 1.37173181e-1, 1.29936675e-1, 1.22838921e-1, 1.15833617e-1, 1.08922157e-1,
      1.02118262e-1, 9.54374347e-2, 8.88965361e-2, 8.25343122e-2, 7.64613376e-2, 7.07768506e-2,
      6.55491245e-2, 6.08209513e-2, 5.65919008e-2, 5.27658978e-2, 4.92440039e-2, 4.59470203e-2,
      4.28122756e-2, 3.98013706e-2, 3.69227010e-2, 3.41908878e-2, 3.16158257e-2, 2.92032829e-2,
      2.69535703e-2, 2.48575196e-2, 2.29046814e-2, 2.10854954e-2, 1.93911943e-2, 1.78145088e-2,
      1.63514395e-2, 1.49980399e-2, 1.37496906e-2, 1.26012838e-2, 1.15465960e-2, 1.05765907e-2,
      9.68267908e-3, 8.85745717e-3, 8.09453368e-3, 7.38889974e-3, 6.73798703e-3, 6.13946693e-3,
      5.59073481e-3, 5.08900051e-3, 4.63119668e-3, 4.21365825e-3, 3.83286115e-3, 3.48559480e-3,
      3.16893322e-3, 2.88020697e-3, 2.61697507e-3, 2.37701167e-3, 2.15828823e-3, 1.95895533e-3,
      1.77737071e-3, 1.61218479e-3, 1.46214596e-3, 1.32605500e-3, 1.20277126e-3, 1.09118319e-3,
      9.90135012e-4, 8.98563819e-4, 8.15524167e-4, 7.40173919e-4, 6.71772381e-4, 6.09693161e-4,
      5.53372350e-4, 5.02292483e-4, 4.55979308e-4, 4.13996781e-4, 3.75939804e-4, 3.41439153e-4,
      3.10160251e-4, 2.81799951e-4, 2.56083822e-4, 2.32763970e-4, 2.11615771e-4, 1.92435527e-4,
      1.75038547e-4, 1.59257316e-4, 1.44939678e-4, 1.31947783e-4, 1.20156843e-4, 1.09453897e-4,
      9.97366706e-5, 9.09125196e-5, 8.28975884e-5, 7.56159926e-5, 6.89990581e-5, 6.29846143e-5,
      5.75163304e-5, 5.25432225e-5, 4.80191572e-5, 4.39023849e-5, 4.01551079e-5, 3.67430788e-5,
      3.36352927e-5, 3.08036839e-5, 2.82228417e-5, 2.58697494e-5, 2.37235394e-5, 2.17653077e-5,
      1.99779303e-5, 1.83458912e-5, 1.68551252e-5, 1.54928703e-5, 1.42475548e-5, 1.31086855e-5,
      1.20667447e-5, 1.11130939e-5, 1.02398848e-5, 9.43999216e-6, 8.70694726e-6, 8.03487566e-6,
      7.41844055e-6, 6.85279186e-6, 6.33352032e-6, 5.85661599e-6, 5.41843079e-6))
  cone_s = np.concatenate((cone_s, np.zeros(215, dtype=cone_s.dtype)))

  # Calculate matching functions
  match_x = 1.94735469 * cone_l - 1.41445123 * cone_m + 0.36476327 * cone_s
  match_y = 0.68990272 * cone_l + 0.34832189 * cone_m
  match_z = 1.93485343 * cone_s

  # Integrate intensities against matching functions
  x = np.zeros_like(intensities[0])
  y = np.zeros_like(intensities[0])
  z = np.zeros_like(intensities[0])
  for n, match_lambda in enumerate(match_lambdas):

    # Skip contribution if wavelength is outside ray-tracing results
    if match_lambda < wavelengths[0] or match_lambda > wavelengths[-1]:
      continue

    # Calculate interpolation of intensities
    ind = np.where(wavelengths <= match_lambda)[0][-1]
    if ind == len(wavelengths) - 1:
      ind -= 1
    frac = (match_lambda - wavelengths[ind]) / (wavelengths[ind+1] - wavelengths[ind])
    intensity = ((1.0 - frac) * wavelengths[ind] * intensities[ind] \
        + frac * wavelengths[ind+1] * intensities[ind+1]) / match_lambda

    # Add contribution to integral
    weight = 1.0
    if n == 0 or n == len(match_lambdas) - 1:
      weight = 0.5
    x += weight * match_x[n] * intensity
    y += weight * match_y[n] * intensity
    z += weight * match_z[n] * intensity

  # Combine results
  return np.array((x, y, z)) / np.nanmax(y)

# Function for converting from XYZ to sRGB1
def xyz_to_rgb(vals):
  x = vals[0]
  y = vals[1]
  z = vals[2]
  lr = np.clip(3.2406 * x - 1.5372 * y - 0.4986 * z, 0.0, 1.0)
  lg = np.clip(-0.9689 * x + 1.8758 * y + 0.0415 * z, 0.0, 1.0)
  lb = np.clip(0.0557 * x - 0.2040 * y + 1.0570 * z, 0.0, 1.0)
  sr = np.where(lr <= 0.0031308, 12.29 * lr, 1.055 * lr ** (1.0 / 2.4) - 0.055)
  sg = np.where(lg <= 0.0031308, 12.29 * lg, 1.055 * lg ** (1.0 / 2.4) - 0.055)
  sb = np.where(lb <= 0.0031308, 12.29 * lb, 1.055 * lb ** (1.0 / 2.4) - 0.055)
  return np.array((sr, sg, sb))

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('filename_data', help='name of file containing raw image data')
  parser.add_argument('filename_plot', help='name of image file to write')
  parser.add_argument('-f', '--central_frequency', type=float,
      help='frequency in Hz to be shifted to center of visual range')
  parser.add_argument('--central_wavelength', type=float, default=550.0,
      help='wavelength in nm to which central_frequency should be shifted')
  parser.add_argument('-d', '--distance', type=float, help='distance to black hole in parsecs')
  parser.add_argument('-a', '--axes', choices=('pixel','rg','cm','muas'), help='axes labels')
  parser.add_argument('-l', '--max_level', type=int, help='maximum adaptive level to plot')
  parser.add_argument('--notex', action='store_true',
      help='flag indicating external tex distribution should not be used')
  args = parser.parse_args()
  main(**vars(args))
