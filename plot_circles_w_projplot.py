import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord

def plot_circle(rot, r, radec=False, fmt="--k"):
  """
  v1.1 mft 12 07 23

  Plots a non-filled circle on top of a HEALPix map

  rot: center of the circle

  r: radius of the circle (in degrees)

  radec: assumes rot is in J2000 coordinates (default: False)
  """

  theta = np.linspace(0, 2*np.pi, 1000)
  lon = rot + r*np.cos(theta)
  lat = rot + r*np.sin(theta)
  hp.projplot(lon, lat, fmt=fmt, lonlat=True)
  return lon, lat
