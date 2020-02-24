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
        her its the simple multiplication without any optimisation
        
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




def random_elliptique(p):
    #while we didnt find a curve
    while(True):
        m=hashlib.sha1()
        seed=secrets.token_bytes(20) #arbitary bit string of 160bits
        m.update(seed)
        hash1=m.digest() #sha1 of seedE
        log2 = math.log(p, 2.0)
        t=math.ceil(log2)
        s=math.floor((t-1)/160)
        h=t-(160*s)
        c0=hash1[0:h] #the h rightmost bit of hash1
        w0=zero_octet_pod_fort(c0) # c0 with the leftmost bit to 0
        W=[]
        W.append(w0)
        #Forming the w it's the concatenation of the Wi 
        #with Wi=sha1((seed+1)mod 2^160)
        for i in range (1,s+1):
            m=hashlib.sha1()
            seedE=((seed+i )%( 2 **160))
            m.update(seedE)
            hash2=m.digest()
            W.append(hash2)
        Wbis=concatenation(W)
        r=0
        #c
        for i in range(1,t+1):
            r=r+(Wbis[len(Wbis)-i]*(2**(t-i)))
        while(True):
            a=random.randint(0,p-1)
            b=random.randint(0,p-1)
            if(((r*(b**2)) % p) == ((a**3) %p)):
                break
        #if we pass this test it mean we have a good curve        
        if((4* (a**3))+(27*(b**2)) % p != 0):
            break
    
    return seed,a,b

def verfy_curve(p,seed,a,b):
    """
    Fonction return true if we generate a good curve or not
    
    """
    m=hashlib.sha1()
    m.update(seed)
    hash1=m.digest()
    log2 = math.log(p, 2.0)
    t=math.ceil(log2)
    s=math.floor((t-1)/160)
    h=t-(160*s)
    c0=hash1[0:h]
    w0=zero_octet_pod_fort(c0)
    W=[]
    W.append(w0)
    for i in range (1,s+1):
        m=hashlib.sha1()
        seedE=((seed+i )%( 2 **160))
        m.update(seedE)
        hash2=m.digest()
        W.append(hash2)
    Wbis=concatenation(W)
    r=0
    for i in range(1,t+1):
        r=r+(Wbis[len(Wbis)-i]*(2**(t-i)))
    if(((r*(b**2)) % p) == ((a**3) %p)):
        return True
    else : 
        return False



def random_point(a,b,p):
    while(True):
        m=hashlib.sha1()
        seed=secrets.token_bytes(20)
        m.update(seed)
        hash1=m.digest()
        log2 = math.log(p, 2.0)
        t=math.ceil(log2)
        s=math.floor((t-1)/160)
        h=t-(160*s)
        c0=hash1[0:h]
        x0=zero_octet_pod_fort(c0)
        X=[]
        X.append(x0)
        for i in range (1,s+1):
            m=hashlib.sha1()
            seedP=((seed+i )%( 2 **160))
            m.update(seedP)
            hash2=m.digest()
            X.append(hash2)
        Xbis=concatenation(X)
        xu=int.from_bytes(Xbis,"little")
        temp=((xu**3)+(a*xu)+b)%p
        if(solution(temp,p)==True):
            break
    for i in range(p):
        if((i**2)%p == temp):
            break
        if(i == p-1):
            print('  Error !!!  ')
    yu=i
    u=Point(EllipticCurve(a,b,p),xu,yu)   
    P=u.__mul__(h)
    return seed,yu,P  

def verify_point(a,b,p,seed,yu,P,ordreP):
    m=hashlib.sha1()
    m.update(seed)
    hash1=m.digest()
    log2 = math.log(p, 2.0)
    t=math.ceil(log2)
    s=math.floor((t-1)/160)
    h=t-(160*s)
    c0=hash1[0:h]
    x0=zero_octet_pod_fort(c0)
    X=[]
    X.append(x0)
    for i in range (1,s+1):
        m=hashlib.sha1()
        seedP=((seed+i )%( 2 **160))
        m.update(seedP)
        hash2=m.digest()
        X.append(hash2)
    Xbis=concatenation(X)
    xu=int.from_bytes(Xbis,"little")
    temp=((xu**3)+(a*xu)+b)%p    
    if(((yu**2) % p) != temp):
        return False
    u=Point(EllipticCurve(a,b,p),xu,yu) 
    Pprime=u.__mul__(h)
    if((P!=Pprime) or( not isinstance(P.__mul__(ordreP),Inf))):
        return False
    return True


#methode take b and transform it with the leftmost bit to 0
def zero_octet_pod_fort(b):
    a=b'\x00'
    l=[]
    for i in range(len(b)):
        if i==len(b)-1:
            l.append(a[0] &b[i])
        else:
            l.append(b[i])
    return bytes(l)

def concatenation(l):
    w=[]
    for wi in l : 
        for i in range(len(wi)):
            w.append(wi[i])
    return bytes(w)


#verify if n is a perfect power in Fp 
def solution(n,p):
    if((n**((p-1)/2))%p==1):
        return True
    else :
        return False


#order of a Point P in Fp
def order(P,p):
    if(isinstance(P,Inf)):
        return 1
    else :
        i=2
    while(True):
        h=P.__mul__(i)
        if(isinstance(h,Inf)):
            return i
        i=i+1
        



x,y,z=random_elliptique(23)  
crouve=EllipticCurve(y,z,23)  
print('seedE value')
print(x)
print(crouve)
print(verfy_curve(23,x,y,z))
    
sed,ye,point=random_point(1,1,23)
print('on a generer ce point man')
print(point.__getitem__(0))
print(point.__getitem__(1))

print(verify_point(1,1,23,sed,ye,point,order(point,23)))


