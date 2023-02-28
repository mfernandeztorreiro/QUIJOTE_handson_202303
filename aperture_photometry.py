import numpy as np
import healpy as hp
from .miscellaneous_functions import Tthermo_factor
Jy = 1e-26

def aperture_photometry(map, freq, cal, pix_object, pix_BG, units="K"):
  """
  v1.1 mft 17 02 23

  Computes aperture photometry estimate for an object in a HEALPix map

  map: input map to perform aperture photometry from

  freq: frequency of the map (in GHz)

  cal: uncertainty calibration factor (0 < cal < 1)

  pix_object: object aperture pixel selection (HEALPix RING order) (must be an array)

  pix_BG: background aperture pixel selection (HEALPix RING order) (must be an array)

  units: must be "K" (temperature; in thermodynamic units) or "MJy/sr" (intensity)
  
  RETURNS

  Sf: extracted flux density (in Jy)
  Serror: its estimated uncertainty (already taking into account 
    calibration uncertainty)
  """
  d2r = np.pi/180
  omegapix = 4 * np.pi / len(map)
  sa = omegapix * len(pix_object)
  nbeams = hp.npix2nside(len(map)) / 64

  if units=="MJy/sr":
    int_BG = np.median(1e6*map[list(pix_BG)])
    int_object = np.mean(1e6*map[list(pix_object)])

    Sf = (int_object - int_BG) * sa
    Serror = np.std(1e6*map[list(pix_BG)]) * sa * np.sqrt(
      (np.pi/2)*nbeams/len(pix_BG) + nbeams/len(pix_object))

  elif units=="K":
    T_BG = np.median(map[list(pix_BG)])
    T_object = np.mean(map[list(pix_object)])
    Tc = T_object - T_BG

    Sf = Tthermo_factor(freq) * Tc / Jy * sa
    Serror = np.std(map[list(pix_BG)]) * Tthermo_factor(freq) * np.sqrt(
      (np.pi/2)*nbeams/len(pix_BG) + nbeams/len(pix_object)) / Jy * sa

  Serror = np.sqrt(Serror**2 + (cal * Sf)**2)

  return Sf, Serror
