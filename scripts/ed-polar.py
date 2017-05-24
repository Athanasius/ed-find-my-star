#!/usr/bin/env python3
#
#  Take a J2000 RA/Dec position and convert into RA/Dec in the ED galaxy
#  frame of reference.

import numpy as np
import math as m

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

# FK5 RA is in hours, minutes, seconds, so first convert to decimal hours
# then divide by 24.0 to get range 0.0 to 1.0, then multiply by 360.0 for
# degrees.
inRA = (14.0 + (39.0 / 60.0) + (36.204 / 3600.0)) / 24.0 * 360.0
# FK5 Dec is already in degrees, minutes, seconds
inDec = -60 - (50.0 * 60.0) - (8.23 / 3600.0)
inParallax = 0.742
lightYearsPerParsec = 3.261563777

def polar2cartesian(ra, dec, length):
  x = length * m.sin(m.radians(dec)) * m.cos(m.radians(ra))
  y = length * m.sin(m.radians(dec)) * m.sin(m.radians(ra))
  z = length * m.cos(m.radians(dec))

  return (x, y, z)

def main():
  distance = 1 / inParallax * lightYearsPerParsec

  (x, y, z) = polar2cartesian(inRA, inDec, distance)

  print("x: %7.5f" % x)
  print("y: %7.5f" % y)
  print("z: %7.5f" % z)

if __name__ == '__main__':
  main()
