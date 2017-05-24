#!/usr/bin/env python3
#
#  Take a J2000 RA/Dec position and convert into RA/Dec in the ED galaxy
#  frame of reference.

import math as m
import numpy as np

# Alpha Centauri, according to simbad:
#
#     FK5 coord. (ep=J2000 eq=2000) :	14 39 36.204 -60 50 08.23 [ 1500 1250 179 ]
#     Gal coord. (ep=J2000) : 315.7330 -00.6809 [ 1500 1250 179 ]
#	    Proper motions mas/yr :	-3608 686 [30 25 179] C ~
#     Parallaxes (mas): 742 [~] D ~
#
# And in ED:
#
# |      name      |    x    |    y     |    z    |
# | Alpha Centauri | 3.03125 | -0.09375 | 3.15625 |

edX = 3.03125
edY = -0.09375
edZ = 3.15625
# FK5 RA is in hours, minutes, seconds, so first convert to decimal hours
# then divide by 24.0 to get range 0.0 to 1.0, then multiply by 360.0 for
# degrees.
inRA = (14.0 + (39.0 / 60.0) + (36.204 / 3600.0)) / 24.0 * 360.0
# FK5 Dec is already in degrees, minutes, seconds, -1.0 here is simply
# because it's negative
inDec = -1.0 * (60.0 + (50.0 / 60.0) + (8.23 / 3600.0))
inParallax = 0.742
lightYearsPerParsec = 3.261563777
inR = 1 / inParallax * lightYearsPerParsec

def sphericalToCartesian(r, ra, dec):
  # We have 'Declination' which is with respect to the equator, being +/-
  # 90 degrees.  The following equation is based on traditional spherical
  # co-ordinates where this is the angle down from the zenith (north pole)
  # So  90 ->  90 - 90 =    0 -> * -1.0 =  0
  #     80 ->  80 - 90 =  -10 -> * -1.0 = 10
  #      0 ->   0 - 90 =  -90 -> * -1.0 = 90
  #    -90 -> -90 - 90 = -180 -> * -1.0 = 180
  dec = (dec - 90.0) * -1.0
  x = r * m.sin(m.radians(dec)) * m.cos(m.radians(ra))
  y = r * m.sin(m.radians(dec)) * m.sin(m.radians(ra))
  z = r * m.cos(m.radians(dec))

  return (x, y, z)

def cartesianToSpherical(x, y, z):
  r = m.sqrt(x ** 2 + y ** 2 + z ** 2)
  print("cartesianToSpherical(): r = %8.5f" % r)
  print("cartesianToSpherical(): z = %8.5f" % z)
  dec = m.degrees(m.acos(z / r))
  ra = m.degrees(m.atan2(y, x))
  print("cartesianToSpherical(): dec (az = zero) = %8.5f" % dec)
  dec = (90.0 - dec)

  return(r, ra, dec)

###########################################################################
# From <https://stackoverflow.com/a/17838396>
###########################################################################
def project_points(x, y, z, a, b, c):
    """
    Projects the points with coordinates x, y, z onto the plane
    defined by a*x + b*y + c*z = 1
    """
    vector_norm = a*a + b*b + c*c
    normal_vector = np.array([a, b, c]) / np.sqrt(vector_norm)
    point_in_plane = np.array([a, b, c]) / vector_norm

    points = np.column_stack((x, y, z))
    points_from_point_in_plane = points - point_in_plane
    proj_onto_normal_vector = np.dot(points_from_point_in_plane,
                                     normal_vector)
    proj_onto_plane = (points_from_point_in_plane -
                       proj_onto_normal_vector[:, None]*normal_vector)

    return point_in_plane + proj_onto_plane
###########################################################################

###########################################################################
# From <https://stackoverflow.com/a/13849249>
###########################################################################
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
###########################################################################

def main():
  edDistance = m.sqrt(edX ** 2 + edY ** 2 + edZ ** 2)
  distance = 1 / inParallax * lightYearsPerParsec

  # Transform FK5 spherical into ED galaxy spherical
  # In 2D this would be easy, just find the angle of the 2D vector to the +ve
  # x vector (1, 0) for both, and then you know what angle to subtract/add to
  # transform from one co-ordinate system to the other.
  #
  # In 3D we can start by doing the same. Project both 3D vectors onto the
  # x/y plane, then find the angle between those and the (1, 0) vector.  This
  # tells us the necessary change in RA between the two systems.
  #
  # THEN just look at the difference in Dec to know how to transform that!
  #
  # So first convert both into spherical:
  #  NB: ED co-ordinate system, viewed from Sol towards Sag A* is
  #          x = right/left, y = up/down, z = forwards/backwards
  #      Maths will have:
  #          x = right/left, y = forwards/backwards, z = up/down
  #      So flip edY and edZ here
  (edR, edRA, edDec) = cartesianToSpherical(edX, edZ, edY)
  print("ED       R: %8.5f   RA: %8.5f  Dec: %8.5f" % (edR, edRA, edDec))
  (check_edX, check_edY, check_edZ) = sphericalToCartesian(edR, edRA, edDec)
  print("ED       x: %8.5f    y: %8.5f    z: %8.5f" % (edX, edY, edZ))
  # Again ED Y / Z flipped to match back up
  print("EDcheck  x: %8.5f    y: %8.5f    z: %8.5f" % (check_edX, check_edZ, check_edY))
  # Astronomical is already in spherical
  # Convert spherical to cartesian (edDistance, not distance, to ensure it
  # matches given the in-game 1/32 ly grid)
  print("IN       R: %8.5f   RA: %8.5f  Dec: %8.5f" % (inR, inRA, inDec))
  (inX, inY, inZ) = sphericalToCartesian(edDistance, inRA, inDec)
  print("IN       x: %8.5f    y: %8.5f    z: %8.5f" % (inX, inY, inZ))
  # Now project both onto x/y plane
  #ed_xy = project_points(edX, edY, edZ, 1, 1, 0)
  #in_xy = project_points(inX, inY, inZ, 1, 1, 0)
  ed_xy = np.array([edX, edZ, 0])
  print(ed_xy)
  in_xy = np.array([inX, inY, 0])
  print(in_xy)
  # Angles between both and (1, 0 , 0)
  ed_xy_angle = m.degrees(angle_between(ed_xy, np.array([1, 0, 0])))
  in_xy_angle = m.degrees(angle_between(in_xy, np.array([1, 0, 0])))
  delta_ra = ed_xy_angle - in_xy_angle
  # And the delta between Decs
  delta_dec = edDec - inDec

  #print(ed_xy)
  #print(in_xy)
  #print("ED -> In RA  angle = %8.5f" % delta_ra)
  #print("edDec = %8.5f" % edDec)
  #print("inDec = %8.5f" % inDec)
  #print("ED -> In Dec angle = %8.5f" % delta_dec)

  outR = inR
  outRA = inRA + delta_ra
  outDec = inDec + delta_dec
  (outX, outY, outZ) = sphericalToCartesian(outR, outRA, outDec)

  print("OUT      x: %8.5f    y: %8.5f    z: %8.5f" % (outX, outY, outZ))
  print("edDistance: %8.5f\tdistance: %8.5f" % (edDistance, distance))

if __name__ == '__main__':
  main()
