"""

"""
from __future__ import annotations # PEP 563: Postponed evaluation of annotations
# since 3.11.4 so no need to use future annotation feature which redefined in python 3.5
# edit: need to use because "__add__(self, other: Point) -> Point:"
from dataclasses import dataclass

#------------------------------------------------------------------------------#
#-----------------------------For wild import----------------------------------#
__all__ = ['Curve', 'Point', 'Generator']

#------------------------------------------------------------------------------#
#-----------------------------Helper functions---------------------------------#

def extended_euclidean_algorithm(a, b):
    """
    Returns (gcd, x, y) s.t. a * x + b * y == gcd
    This function implements the extended Euclidean
    algorithm and runs in O(log b) in the worst case,
    taken from Wikipedia.
    """
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    return old_r, old_s, old_t

def inv(n, p):
    """ returns modular multiplicate inverse m s.t. (n * m) % p == 1 """
    gcd, x, y = extended_euclidean_algorithm(n, p) # pylint: disable=unused-variable
    return x % p


#------------------------------------------------------------------------------#
#-----------------------------Main objects-------------------------------------#

@dataclass
class Curve:
    """
    Elliptic Curve over the field of integers modulo a prime.
    Points on the curve satisfy y^2 = x^3 + a*x + b (mod p).

    Z_p là một TẬP HỢP, khi p là số nguyên tố thì có thể coi nó là một TRƯỜNG https://forum.mathscope.org/archive/index.php/t-11835.html
    """
    p: int # ta nói secp256k1 có đặc trưng p, được định trong trường Z_p
    a: int
    b: int


@dataclass
class Point:
    """ Số nguyên tọa độ (x,y) trên đường cong """
    curve: Curve
    x: int
    y: int

    def __add__(self, other: Point) -> Point:
        # handle special case of P + 0 = 0 + P = 0
        if self == INF:
            return other
        if other == INF:
            return self
        # handle special case of P + (-P) = 0
        if self.x == other.x and self.y != other.y:
            return INF
        # compute the "slope"
        if self.x == other.x: # (self.y = other.y is guaranteed too per above check)
            m = (3 * self.x**2 + self.curve.a) * inv(2 * self.y, self.curve.p)
        else:
            m = (self.y - other.y) * inv(self.x - other.x, self.curve.p)
        # compute the new point
        rx = (m**2 - self.x - other.x) % self.curve.p
        ry = (-(m*(rx - self.x) + self.y)) % self.curve.p
        return Point(self.curve, rx, ry)
    
    def __rmul__(self, other: Point) -> Point:
        assert isinstance(k, int) and k >= 0
        result = INF
        append = self
        while k:
            if k & 1:
                result += append
            append += append
            k >>= 1
        return result

@dataclass
class Generator:
    """
    A Generator over a curve: an initial point and pre-declared order.
    """
    G: Point # starting point on the curve
    n: int # the order 0*G = n*G = INF (????)


INF = Point(None, None, None)