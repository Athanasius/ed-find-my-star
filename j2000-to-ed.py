#!/usr/bin/env python3
#
#
#

import sys
import math
import numpy as np

np.set_printoptions(12)
# This is the 'Nj' array from https://arxiv.org/pdf/1010.3773.pdf
nj = np.array([[-0.054875539390, -0.873437104725, -0.483834991775],
               [ 0.494109453633, -0.444829594298,  0.746982248696],
               [-0.867666135681, -0.198076389622,  0.455983794523] ])
# This is the 'Nb' array from https://arxiv.org/pdf/1010.3773.pdf
nb = np.array([[-0.066988739410, -0.872755765850, -0.483538914637],
               [ 0.492728466081, -0.450346958020,  0.744584633279],
               [-0.867600811149, -0.188374601732,  0.460199784785] ])
# And this is me trying to get one to fit in-game data
ne = np.array([[-0.066988739410, -0.872755765850, -0.483538914637],
               [ 0.492728466081, -0.450346958020,  0.744584633279],
               [-0.86747076, -0.17863258,  0.4643112 ] ])
#               [-0.867666135681, -0.17867281384199998,  0.464415775047] ])

# *) Convert input RA/Dec to decimal degrees
# Alpha Centauri:
#   http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=alpha+centauri
# Using FK5 (ep=J2000 eq=2000)
#   ra = 14 39 36.204
#   dec = -60 50 08.23
#   Parallax = 742 (milli-arcseconds)

###########################################################################
# Right Ascension from first three command-line arguments
###########################################################################
ra = {'hour': float(sys.argv[1]), 'minute': float(sys.argv[2]), 'second': float(sys.argv[3])}
#ra = {'hour': 14, 'minute': 39, 'second': 36.204}
ra['degrees'] = (ra['hour'] + ra['minute']/60.0 + ra['second']/60.0/60.0) / 24.0 * 360.0
###########################################################################

###########################################################################
# Declination from 4th through 6th command-line arguments
###########################################################################
if sys.argv[4][0] == '+':
  dec_degrees = sys.argv[4][1:]
else:
  dec_degrees = sys.argv[4]
dec = {'degree': float(dec_degrees), 'minute': float(sys.argv[5]), 'second': float(sys.argv[6])}
#dec = {'degree': -60, 'minute': 50, 'second': 8.23}
# Declination has a sign, but is only signified on the 'degrees' field, so
# if it's not 0 we use Dec_Degrees / abs(Dec_Degrees) to store the sign
# store abs(Dec_Degrees) back in it, and then multiply by the sign afterwards
if dec['degree'] != 0.0:
  dec['sign'] = dec['degree'] / abs(dec['degree'])
else:
  dec['sign'] = 1.0
dec['degree'] = abs(dec['degree'])
#dec = {'degree': -60, 'minute': 50, 'second': 8.23}
dec['degrees'] = dec['sign'] * ( dec['degree'] + dec['minute']/60.0 + dec['second']/60.0/60.0 )
###########################################################################

###########################################################################
# Parallax (milli-arcseconds) from the 7th command-line argument
###########################################################################
#parallax = 742.0
parallax = float(sys.argv[7])
###########################################################################

### HIP 71181
#ra = (14.0 + 33.0/60.0 + 28.868/60.0/60.0) / 24.0 * 360.0
#dec = 52.0 + 54.0/60.0 + 31.65/60.0/60.0
#parallax = 75.65
# Convert to radians
ra['rads'] = math.radians(ra['degrees'])
dec['rads'] = math.radians(dec['degrees'])
print("RA :", ra)
print("Dec: ", dec)
# *) Xeu = cos(ra) * cos(dec)
#    Yeu = sin(ra) * cos(dec)
#    Zeu = sin(dec)
eu = np.array([math.cos(ra['rads']) * math.cos(dec['rads']),
               math.sin(ra['rads']) * math.cos(dec['rads']),
               math.sin(dec['rads']) ])
print("EU        ", eu)

# *) Multiply our vector by the Nj array
print("Using nb")
gu = np.matmul(nb, eu)
print("GU        ", gu)

# *) We have a unit vector, so multiply up by the distance.
#    parsecs = 1/<parallax in seconds>
#    lightyears = parsecs * 3.26

parsecs = 1 / (parallax / 1000.0)
print("Distance (Parsecs): %f" % (parsecs))
lightyears = parsecs * 3.26
print("Distance (Lightyears): %f" % (lightyears))

guDistant = gu * lightyears
print("guDistant ", guDistant)
ed = np.array([-1.0 * guDistant[1], guDistant[2], guDistant[0]])

print("ED        ", ed)

# Normalise to the 1/32 grid ED effectively uses
ed[0] = int(ed[0] * 32.0) / 32.0
ed[1] = int(ed[1] * 32.0) / 32.0
ed[2] = int(ed[2] * 32.0) / 32.0
print("ED(normal)", ed)

#      name      |    x    |    y     |    z    | argv
#----------------+---------+----------+---------+-------------------------
# Alpha Centauri | 3.03125 | -0.09375 | 3.15625 | 14 39 36.204 -60 50 08.23 742.0
# ED(normal) [ 3.03125 -0.09375  3.15625]
#      Ross 128  | 5.53125 |  9.4375  |   0.125 | 11 47 44.397 00 48 16.40 295.80
# 
#      Luhman 16 |  6.3125 | 0.59375  | 1.71875 | not in simbad?
# Barnard's Star | -3.03125|    1.375 |   4.9375| 17 57 48.498 +04 41 36.21 548.31
#
