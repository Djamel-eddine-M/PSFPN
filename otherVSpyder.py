#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:31:03 2020

@author: 3673868
"""

import secrets
import math
import hashlib
import random
import time


class EllipticCurve(object):
    """Represents a single elliptic curve defined over a finite field.

    See here:
        http://en.wikipedia.org/wiki/Elliptic_curve
        

    p must be prime, since we use the modular inverse to compute point
    addition.

    """
    def __init__(self, a, b, p):
        self.a = a%p
        self.b = b%p
        self.p = p

    def __eq__(self, C):
        return (self.a, self.b) == (C.a, C.b)

    def has_point(self, x, y):
        return (y ** 2) % self.p == (x ** 3 + self.a * x + self.b) % self.p

    def __str__(self):
        return 'y^2 = x^3 + {}x + {}'.format(self.a, self.b)



class Point(object):
    """A point on a specific curve. all points are associated with a curve
        so we have to test if this point are in the curve 
        and we have another class for the infinite point
    """
    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x % curve.p
        self.y = y % curve.p
        #we test if the point are on the curve
        if not self.curve.has_point(x, y):
            raise ValueError('{} is not on curve {}'.format(self, self.curve))

    def __str__(self):
        return '({}, {})'.format(self.x, self.y)

    def __getitem__(self, index):
        return [self.x, self.y][index]

    def __eq__(self, Q):
        return (self.curve, self.x, self.y) == (Q.curve, Q.x, Q.y)

    def __neg__(self):
        return Point(self.curve, self.x, -self.y)
    
    def __add__(self, Q):
        """Add two points together.

        for that we need to take care of special cases:
         * Q is the infinity point (0)
         * P == Q
         * The line crossing P and Q is vertical.

        """
        assert self.curve == Q.curve
        # 0 + P = P
        if isinstance(Q, Inf):
            return self

        xp, yp, xq, yq = self.x, self.y, Q.x, Q.y
        m = None

        # P == Q
        if self == Q:
            if self.y == 0:
                R = Inf(self.curve)
            else:
                m = ((3 * xp * xp + self.curve.a) * mod_inverse(2 * yp, self.curve.p)) % self.curve.p

        # Vertical line
        elif xp == xq:
            R = Inf(self.curve)

        # Common case
        else:
            m = ((yq - yp) * mod_inverse(xq - xp, self.curve.p)) % self.curve.p
        #if it's not the point at the infinite
        if m is not None:
            xr = (m ** 2 - xp - xq) % self.curve.p
            yr = (m * (xp - xr) - yp) % self.curve.p
            R = Point(self.curve, xr, yr)


        return R
    
    
    def __mul__(self, n):
        
        """ Multiplication of a point by a integer
           
        """
        assert isinstance(n, int)
        assert n > 0

        n = n % self.curve.p
    
        if n == 0:
            return Inf(self.curve)
    
        else:
            Q = self
            R = Inf(self.curve)
    
            i = 1
            while i <= n:
                if n & i == i:
                    R = R + Q
    
                Q = Q + Q
    
                i = i << 1
    
        return R

    def __rmul__(self, n):
        return self * n
    
    
class Inf(Point):
    """The custom infinity point."""
    def __init__(self, curve):
        self.curve = curve

    def __eq__(self, Q):
        return isinstance(Q, Inf)

    def __neg__(self):
        """-0 = 0"""
        return self
    def __add__(self, Q):
        """P + 0 = P"""
        return Q
    
    def __str__(self):
        return 'Inf'

def mod_inverse(a, n):
    """Return the inverse of a mod n.

    n must be prime.

    >>> mod_inverse(42, 2017)    


    """
    b = n
    if abs(b) == 0:
        return (1, 0, a)

    x1, x2, y1, y2 = 0, 1, 1, 0
    while abs(b) > 0:
        q, r = divmod(a, b)
        x = x2 - q * x1
        y = y2 - q * y1
        a, b, x2, x1, y2, y1 = b, r, x1, x, y1, y

    return x2 % n



def random_elliptiqueV2(p):
    #while we