import numpy as np
import healpy as hp

def apertures_definition(nside, rot=[184.55896,-5.78467], r1=1.5, r2=None,
  theta=[0,360], rot_BG=None, r3=2, inclusive=True):
  """
  v1.1 mft 17 02 23
  Defines apertures required to perform aperture photometry studies. 
  Default values for CRAB / TAU A

  nside: HEALPix factor, defines the number of equal area pixels over
    which the sky sphere is distributed

  rot: [theta,phi] position of the desired object to study in galactic
    coordiantes (default: crab position) (iterable)

  rot_BG: center of the BG region (default: same as rot) (interable)

  r1: radius of the primary (object) aperture, in degrees

  r2: inner radius of the BG region (centered on rot_BG) (default:
    equal to r1), in degrees

  r3: outer raidus of the BG region (centered on rot_BG), in degrees

  theta: defines sections of the previous ring as the BG. E.g. theta=[120,170]
    -> only the section between 120 and 170 degrees (with 0 being the South and
    180 the North, clockwise) is considered to estimate the BG level.
    (iterable: length must be even)
  
  inclusive: value for "inclusive" parameter in hp.query_disc or hp.query_polygon (default: True)

  RETURNS

  listpix_object: 1D array with the indeces for the pixels inside the primary
    aperture

  ring_selection: same but for the BG region
  """

  rot_BG = rot if rot_BG is None else rot_BG
  r2 = r1 if r2 is None else r2

  d2r = np.pi / 180
  x, y, z= hp.ang2vec((90-rot[1])*d2r, rot[0]*d2r)
  x_BG, y_BG, z_BG= hp.ang2vec((90-rot_BG[1])*d2r, rot_BG[0]*d2r)
  listpix_object = hp.query_disc(nside=nside, vec=[x,y,z], radius=r1*d2r,
    inclusive=inclusive)
  listpix_r2 =hp.query_disc(nside, vec=[x_BG,y_BG,z_BG], radius=r2*d2r,
    inclusive=inclusive)
  listpix_r3 =hp.query_disc(nside, vec=[x_BG,y_BG,z_BG], radius=r3*d2r,
    inclusive=inclusive)
  ring = set(listpix_r3) - set(listpix_r2)

  if theta==[0,360]:
    return listpix_object, np.array([list(ring)])[0]

  ring_selection = []
  for i in ring:
    pix_pos = np.array(hp.pix2ang(nside, i)) - np.array([(90-rot_BG[1])*d2r,
      rot_BG[0]*d2r])
    c = pix_pos[0] / np.sqrt(np.sum(pix_pos**2))
    s = pix_pos[1] / np.sqrt(np.sum(pix_pos**2))
    if s<0:
      alpha = 2*np.pi - np.arccos(c)
    else:
      alpha = np.arccos(c)

    for j in range(int(len(theta)/2)):
      if alpha>theta[2*j]*d2r and alpha<theta[2*j+1]*d2r:
        ring_selection.append(i)

  return listpix_object, np.array([list(ring_selection)])[0]