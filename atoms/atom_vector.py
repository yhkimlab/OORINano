#
# NanoCore3
# Last revision         : 2019. 02. 12
# modularized by J. Park: 2021. 10. 30
#

from __future__ import print_function
#from .atomic_data import atomic_weight, atomic_symbol, covalent_radii 
from .atomic_data import *
from math import sqrt, pi, sin, cos, asin, acos
import numpy as np
import copy     # added by SH
import os
import sys

## Vector #

class Vector(object):
    """
    class Vector to replace Scientific.Geometry.Vector
    """
    
    __slots__ = ['x1', 'x2', 'x3']
    
    def __init__(self, x1, x2=None, x3=None):
        if x2 == None and x3 == None:
            self.x1 = float(x1[0]); self.x2 = float(x1[1]); self.x3 = float(x1[2])
        else:
            self.x1 = float(x1); self.x2 = float(x2); self.x3 = float(x3)

    def __repr__(self):
        return "Vector(%10.6f, %10.6f, %10.6f)" % \
	(self.x1, self.x2, self.x3)

    def __getitem__(self, i):
        return [self.x1, self.x2, self.x3][i]

    def __mul__(self, other):
        """
        Return a cross products vector
        """
        if isinstance(other, Vector):
            a1, a2, a3 = [self.x1, self.x2, self.x3]
            b1, b2, b3 = [other.x1, other.x2, other.x3]
            c1 =   a2*b3 - a3*b2
            c2 = -(a1*b3 - a3*b1)
            c3 =   a1*b2 - a2*b1
            return Vector(c1, c2, c3)
        else:
            return Vector(other*self.x1, other*self.x2, other*self.x3)

    def __rmul__(self, other):
        """
        Return a cross products vector
        """
        return Vector(other*self.x1, other*self.x2, other*self.x3)

    def cross(self, other): return self*other

    def dot(self, other):
        """
        Return a dot products value
        """
        a1, a2, a3 = [self.x1, self.x2, self.x3]
        b1, b2, b3 = [other.x1, other.x2, other.x3]
        return a1*b1 + a2*b2 + a3*b3

    def x(self): return self.x1
    def y(self): return self.x2
    def z(self): return self.x3

    def __add__(self, other):
        a1, a2, a3 = [self.x1, self.x2, self.x3]
        b1, b2, b3 = [other.x1, other.x2, other.x3]
        return Vector(a1+b1, a2+b2, a3+b3)

    def __sub__(self,other):
        a1, a2, a3 = [self.x1, self.x2, self.x3]
        b1, b2, b3 = [other.x1, other.x2, other.x3]
        return Vector(a1-b1, a2-b2, a3-b3)

    def __div__(self, other):
        other = float(other)
        return Vector(self.x1/other, self.x2/other, self.x3/other)

    def __truediv__(self, other):
        other = float(other)
        return Vector(self.x1/other, self.x2/other, self.x3/other)

    def __abs__(self): return self.length()

    def __len__(self): return 3
    
    def length(self): return sqrt(self.x1**2 + self.x2**2 + self.x3**2)

    def angle(self,other):
        """
        return angle between vector1 and vector2 in radian unit
        >>> angle = vector1.angle(vector2)
        """
        adotb = self.dot(other)
        lena  = self.length(); lenb = other.length()
        costh = adotb / (lena*lenb)
        return acos(costh)
  
    def rotate(self, angle, rot_center, rot_vector):
        """
        rotate a vector using Rodrigues' rotation formula
        >>> angle = 30.0 # in degree
        >>> rot_center = Vector(0,0,0) # center of rotation
        >>> rot_vector = Vector(0,0,1) # rotation axis
        >>> vector.rotate(angle, rot_center, rot_vector)
        """
        rot_v = rot_vector/rot_vector.length()       # rotation axis vector
        obj_v = Vector(self.x(), self.y(), self.z()) # original vector
        ux, uy, uz = rot_v.x(), rot_v.y(), rot_v.z() # rotation axis vector component
        ax, ay, az = obj_v - rot_center              # original vector component centered at (0,0,0)
        ang_rad = angle*pi/180.                      # rotation angle in radian
        cos_a = cos(ang_rad); sin_a = sin(ang_rad)   # cos, sin value of rotation angle
        umat = np.matrix([ax, ay, az]).T             # matrix form of original vector
        # rotation matrix
        Rmat = np.matrix( [[ cos_a + (ux**2)*(1-cos_a),  ux*uy*(1-cos_a) - uz*sin_a, ux*uz*(1-cos_a) + uy*sin_a ], 
                           [ uy*ux*(1-cos_a) + uz*sin_a, cos_a + (uy**2)*(1-cos_a),  uy*uz*(1-cos_a) - ux*sin_a ],
                           [ uz*ux*(1-cos_a) - uy*sin_a, uz*uy*(1-cos_a) + ux*sin_a, cos_a + (uz**2)*(1-cos_a)  ]] )
        # rotated vector componets
        bx, by, bz = Rmat*umat
        # rotated vector vector component centered at the original center
        res_v = Vector(bx, by, bz) + rot_center
        return res_v


