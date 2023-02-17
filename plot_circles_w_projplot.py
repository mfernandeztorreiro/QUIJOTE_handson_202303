import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord

def plot_circle(rot, r, fmt="--k", radec=False, return_coord=False):
  """
  v1.1 mft 12 07 23

  Plots a non-filled circle on top of a HEALPix map

  rot: center of the circle

  r: radius of the circle (in degrees)

  radec: assumes rot is in J2000 coordinates (default: False)
  """

  theta = np.linspace(0, 2*np.pi, 1000)
  lon = rot[0] + r*np.cos(theta)
  lat = rot[1] + r*np.sin(theta)
  hp.projplot(lon, lat, fmt, lonlat=True)
  if return_coord==True:
    return lon, lat
