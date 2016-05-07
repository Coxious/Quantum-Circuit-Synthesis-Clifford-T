# date: May 7, 2016
# brief: EuclideanDomain Z[1/sqrt(2),i] and Z[exp(Pi*I/4)]

from sympy import *
from sympy.polys.rings import *
init_printing(use_unicode=True)

# Ring Number Z[i,1/sqrt(2)] -> Z[w]
class EuclideanDomainNum:
    def __init__(self):
        self.root = exp(pi*I*Rational(1,4))
        self.poly = [0,0,0,0]
        self.de = 0

    def __init__(self, poly, de):
        self.root = exp(pi*I*Rational(1,4))
        self.set(poly, de)

    def set(self, poly, de):
        assert isinstance(poly,list)
        assert len(poly)==4
        assert isinstance(poly[0], int)
        assert isinstance(poly[1], int)
        assert isinstance(poly[2], int)
        assert isinstance(poly[3], int)
        self.poly = poly
        self.de = de

    # Adjust all coefficients after multiplied by sqrt2
    def mul_sqrt2(self, r):
        s = [r[2]-r[0], r[1]+r[3], r[2]+r[0], r[1]-r[3]]
        return s

    def reduce(self, type):
        # Re-calculate its denomExp and poly
        if type == 0:
            # For 2
            isDiv2 = True
            while isDiv2:
                for i in range(0,4):
                    isDiv2 = isDiv2 and (self.poly[i]%2==0)
                if isDiv2:
                    for i in range(0,4):
                        self.poly[i] *= Rational(1,2)
                    self.de -= 2
            # For sqrt2
            isDivSqrt2 = ((self.poly[2]-self.poly[0])%2 == 0) \
                        and ((self.poly[1]+self.poly[3])%2 == 0) \
                        and ((self.poly[2]+self.poly[0])%2 == 0) \
                        and ((self.poly[1]-self.poly[3])%2 == 0)
            if isDivSqrt2:
                poly = [0,0,0,0]
                poly[3] = (self.poly[2]-self.poly[0])*Rational(1,2)
                poly[2] = (self.poly[1]+self.poly[3])*Rational(1,2)
                poly[1] = (self.poly[2]+self.poly[0])*Rational(1,2)
                poly[0] = (self.poly[1]-self.poly[3])*Rational(1,2)
                self.poly = poly
                self.de -= 1
        return self

    def __add__(self, other):
        poly = [0,0,0,0]
        r = self.poly
        s = other.poly
        de = max(self.de, other.de)
        if self.de != other.de:
            if self.de < other.de:
                r = other.poly
                s = self.poly
            times = abs(self.de-other.de)
            while times > 0:
                times -= 1
                s = self.mul_sqrt2(s)
        for i in range(0, 4):
            poly[i] = r[i] + s[i]
        return EuclideanDomainNum(poly, de).reduce(0)

    def __sub__(self, other):
        poly = [0,0,0,0]
        for i in range(0,4):
            other.poly[i] *= -1
        return self+other

    def __mul__(self, other):
        poly = [0,0,0,0]
        r, s = self.poly, other.poly
        de = self.de + other.de
        poly[3] = r[3]*s[3]-r[2]*s[0]-r[1]*s[1]-r[0]*s[2]
        poly[2] = r[2]*s[3]+r[3]*s[2]-r[1]*s[0]-r[0]*s[1]
        poly[1] = r[3]*s[1]+r[2]*s[2]+r[1]*s[3]-r[0]*s[0]
        poly[0] = s[0]*r[3]+s[1]*r[2]+s[2]*r[1]+s[3]*r[0]
        return EuclideanDomainNum(poly, de).reduce(0)

    def __eq__(self, other):
        for i in range(0,4):
            if self.poly[i] != other.poly[i]:
                return False
        return self.de == other.de

    def v(self):
        re = self.poly[3]+sqrt(Rational(1,2))*(self.poly[2]-self.poly[0])
        im = self.poly[1]+sqrt(Rational(1,2))*(self.poly[2]+self.poly[0])
        return ((sqrt(Rational(1,2))**self.de) * (re+I*im)).expand(complex=True)

# Ring number test
a1 = EuclideanDomainNum([2,3,4,5],1)
a2 = EuclideanDomainNum([2,3,4,5],2)
b1 = EuclideanDomainNum([3,4,5,6],1)
b3 = EuclideanDomainNum([3,4,5,6],3)
c0 = EuclideanDomainNum([1,1,1,1],0)
c1 = EuclideanDomainNum([1,1,1,1],1)
print("(1)[2,3,4,5] = ", a1.v())
print("(2)[2,3,4,5] = ", a2.v())
print("(1)[3,4,5,6] = ", b1.v())
print("(3)[3,4,5,6] = ", b3.v())
print("(0)[1,1,1,1] = ", c0.v())
print("(1)[1,1,1,1] = ", c1.v())
print("(1)[2,3,4,5]+(1)[3,4,5,6] = ", (a1+b1).v())
print("(1)[2,3,4,5]+(3)[3,4,5,6] = ", (a1+b3).v())
print("(2)[2,3,4,5]+(1)[3,4,5,6] = ", (a2+b1).v())
print("(1)[2,3,4,5]-(3)[3,4,5,6] = ", (a1-b3).v())
print("(1)[2,3,4,5]+(1)[1,1,1,1] = ", (a1+c1).v())
print("(1)[2,3,4,5]+(1)[1,1,1,1]=(1)[3,4,5,6] : ", a1+c1==b1)
print("(1)[2,3,4,5]+(3)[3,4,5,6]=(0)[1,1,1,1] : ", a1+b3==c0)
print("(2)[2,3,4,5]+(3)[3,4,5,6] = ", (a2*b3).v())