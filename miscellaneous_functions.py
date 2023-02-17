import numpy as np
h = 6.626e-34; kb = 1.38e-23; Tcmb = 2.72548; c=3e8; Jy=1e-26

def x(nu, T=Tcmb):
  """
  v1.1 mft 17 02 23
  Computes x = h * nu / (k_B * T) factor

  nu: frequency to compute (GHz)

  T: temperature (K). Default: CMB temperature (Fixsen+2009)

  RETURNS

  x
  """
  return h * nu * 1e9 / (kb * T)

def Tthermo_factor(nu):
  """
  v1.1 mft 17 02 23
  Computes a(nu) factor from e.g. Génova-Santos+2015, Fernández-Torreiro+2023.
  This factor translates temperature units (in thermodynamic units, K_CMB) to
  intensity units (W / m2 / Hz)

  nu: frequency to compute (GHz)
  """
  return 2 * kb * (nu*1e9)**2 / c**2 * x(nu)**2 * np.exp(x(nu)) / (np.exp(x(nu)) - 1)**2