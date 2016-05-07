# date: May 7, 2016
# brief: Euclidean domain Z[1/sqrt(2),i] and Z[exp(Pi*I/4)]

from sympy import *
from sympy.polys.rings import *
init_printing(use_unicode=True)

# Euclidean Domain Number Z[i,1/sqrt(2)] -> Z[w]
class EuclideanDomainNum:
    def __init__(self):
        self.root = exp(pi*I*Rational(1,4))
        self.poly = [0,0,0,0]
        self.de = 0

    def __init__(self, poly, de):
        self.root = exp(pi*I*Rational(1,4))
        self.set(poly, de)

    # Update its polynomial and denomExp
    def set(self, poly, de):
        assert isinstance(poly,list)
        assert len(poly)==4
        assert isinstance(poly[0], Integer) or isinstance(poly[0], int)
        assert isinstance(poly[1], Integer) or isinstance(poly[1], int)
        assert isinstance(poly[2], Integer) or isinstance(poly[2], int)
        assert isinstance(poly[3], Integer) or isinstance(poly[3], int)
        self.poly = poly
        self.de = de
        self.reduce(0)

    # Adjust all coefficients after multiplied by sqrt2
    def mul_sqrt2(self, r):
        s = [r[1]-r[3], r[2]+r[0], r[1]+r[3], r[2]-r[0]]
        return s

    # Re-calculate its denomExp and poly by case
    def reduce(self, type):
        # Re-calculate its denomExp and poly directly
        if type == 0:
            # For 2
            isDiv2 = True
            while isDiv2 and self.de >= 2:
                for i in range(0,4):
                    isDiv2 = isDiv2 and (self.poly[i]%2==0)
                if isDiv2:
                    for i in range(0,4):
                        self.poly[i] *= Rational(1,2)
                    self.de -= 2
            # For sqrt2
            isDivSqrt2 = ((self.poly[1]-self.poly[3])%2 == 0) \
                     and ((self.poly[2]+self.poly[0])%2 == 0) \
                     and ((self.poly[1]+self.poly[3])%2 == 0) \
                     and ((self.poly[2]-self.poly[0])%2 == 0)
            if isDivSqrt2 and self.de >= 1:
                poly = [0,0,0,0]
                poly[0] = (self.poly[1]-self.poly[3])*Rational(1,2)
                poly[1] = (self.poly[2]+self.poly[0])*Rational(1,2)
                poly[2] = (self.poly[1]+self.poly[3])*Rational(1,2)
                poly[3] = (self.poly[2]-self.poly[0])*Rational(1,2)
                self.poly = poly
                self.de -= 1
        return self

    # Add operation for different donomExp situations
    def __add__(self, other):
        poly = [0,0,0,0]
        r = self.poly
        s = other.poly
        de = max(self.de, other.de)
        # For the smaller one, update its polynomial
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

    # Subtraction operation for different donomExp situations
    def __sub__(self, other):
        poly = [0,0,0,0]
        for i in range(0,4):
            poly[i] = -1 * other.poly[i]
        return self+EuclideanDomainNum(poly, other.de)

    # Multiplication operation for different donomExp situations
    def __mul__(self, other):
        poly = [0,0,0,0]
        r, s = self.poly, other.poly
        de = self.de + other.de
        # Multiplication by its coefficients
        poly[3] = r[3]*s[3]-r[2]*s[0]-r[1]*s[1]-r[0]*s[2]
        poly[2] = r[2]*s[3]+r[3]*s[2]-r[1]*s[0]-r[0]*s[1]
        poly[1] = r[3]*s[1]+r[2]*s[2]+r[1]*s[3]-r[0]*s[0]
        poly[0] = s[0]*r[3]+s[1]*r[2]+s[2]*r[1]+s[3]*r[0]
        return EuclideanDomainNum(poly, de).reduce(0)

    # Equal or not
    def __eq__(self, other):
        for i in range(0,4):
            if self.poly[i] != other.poly[i]:
                return False
        return self.de == other.de

    # Return norm function
    def norm(self, useDenomExp):
        assert isinstance(useDenomExp, bool)
        r = self.poly
        term1 = r[0]**2 + r[1]**2 + r[2]**2 + r[3]**2
        term2 = r[0]*r[1] + r[1]*r[2] + r[2]*r[3] - r[3]*r[0]
        # norm at the Z[exp(Pi*I/4)]
        res = term1**2 - 2*(term2**2)
        if useDenomExp:
            return res*(Rational(1,2)**(2*self.de))
        else:
            return res

    # Return value
    def v(self):
        re = self.poly[3]+sqrt(Rational(1,2))*(self.poly[2]-self.poly[0])
        im = self.poly[1]+sqrt(Rational(1,2))*(self.poly[2]+self.poly[0])
        return ((sqrt(Rational(1,2))**self.de) * (re+I*im)).expand(complex=True)

# Ring number test
a0 = EuclideanDomainNum([2,3,4,5],0)
a1 = EuclideanDomainNum([2,3,4,5],1)
a2 = EuclideanDomainNum([2,3,4,5],2)
b1 = EuclideanDomainNum([3,4,5,6],1)
b3 = EuclideanDomainNum([3,4,5,6],3)
c0 = EuclideanDomainNum([1,1,1,1],0)
c1 = EuclideanDomainNum([1,1,1,1],1)
# Initialization Test
# (1)[2,3,4,5] =  1 + 5*sqrt(2)/2 + 3*sqrt(2)*I/2 + 3*I
# (2)[2,3,4,5] =  sqrt(2)/2 + 5/2 + 3*I/2 + 3*sqrt(2)*I/2
# (1)[3,4,5,6] =  1 + 3*sqrt(2) + 2*sqrt(2)*I + 4*I
# (3)[3,4,5,6] =  1/2 + 3*sqrt(2)/2 + sqrt(2)*I + 2*I
# (0)[1,1,1,1] =  1 + I + sqrt(2)*I
# (1)[1,1,1,1] =  sqrt(2)/2 + sqrt(2)*I/2 + I
print("(0)[2,3,4,5] = ", a1.v(), a1.norm(False))
print("(1)[2,3,4,5] = ", a1.v(), a1.norm(False))
print("(2)[2,3,4,5] = ", a2.v(), a2.norm(False))
print("(1)[3,4,5,6] = ", b1.v(), b1.norm(False))
print("(3)[3,4,5,6] = ", b3.v(), b3.norm(False))
print("(0)[1,1,1,1] = ", c0.v(), c0.norm(False))
print("(1)[1,1,1,1] = ", c1.v(), c1.norm(False))
# Addition Test
# (1)[2,3,4,5]+(1)[3,4,5,6] =  2 + 11*sqrt(2)/2 + 7*sqrt(2)*I/2 + 7*I
# (1)[2,3,4,5]+(3)[3,4,5,6] =  3/2 + 4*sqrt(2) + 5*sqrt(2)*I/2 + 5*I
# (2)[2,3,4,5]+(1)[3,4,5,6] =  7/2 + 7*sqrt(2)/2 + 7*sqrt(2)*I/2 + 11*I/2
# (1)[2,3,4,5]-(3)[3,4,5,6] =  1/2 + sqrt(2) + sqrt(2)*I/2 + I
# (2)[2,3,4,5]-(0)[1,1,1,1] =  sqrt(2)/2 + 3/2 + I/2 + sqrt(2)*I/2
# (1)[2,3,4,5]+(1)[1,1,1,1] =  1 + 3*sqrt(2) + 2*sqrt(2)*I + 4*I
print("(1)[2,3,4,5]+(1)[3,4,5,6] = ", (a1+b1).v())
print("(1)[2,3,4,5]+(3)[3,4,5,6] = ", (a1+b3).v())
print("(2)[2,3,4,5]+(1)[3,4,5,6] = ", (a2+b1).v())
print("(1)[2,3,4,5]-(3)[3,4,5,6] = ", (a1-b3).v())
print("(2)[2,3,4,5]-(0)[1,1,1,1] = ", (a2-c0).v())
print("(1)[2,3,4,5]+(1)[1,1,1,1] = ", (a1+c1).v())
# Equivalent Test
# (1)[2,3,4,5]+(1)[1,1,1,1]=(1)[3,4,5,6] :  True
# (1)[2,3,4,5]+(3)[3,4,5,6]=(0)[1,1,1,1] :  False
print("(1)[2,3,4,5]+(1)[1,1,1,1]=(1)[3,4,5,6] : ", a1+c1 == b1)
print("(1)[2,3,4,5]+(3)[3,4,5,6]=(0)[1,1,1,1] : ", a1+b3 == c0)
# Multiplication Test
# (0)[2,3,4,5]*(0)[1,1,1,1] =  -5*sqrt(2) - 4 + 10*I + 9*sqrt(2)*I
# (1)[2,3,4,5]*(3)[3,4,5,6] =  -13*sqrt(2)/4 - 1 + 13*I + 45*sqrt(2)*I/4
# (2)[2,3,4,5]*(1)[3,4,5,6] =  -13/2 - sqrt(2) + 13*sqrt(2)*I + 45*I/2
# (2)[2,3,4,5]*(3)[3,4,5,6] =  -13/4 - sqrt(2)/2 + 13*sqrt(2)*I/2 + 45*I/4
print("(0)[2,3,4,5]*(0)[1,1,1,1] = ", (a0*c0).v())
print("(1)[2,3,4,5]*(3)[3,4,5,6] = ", (a1*b3).v())
print("(2)[2,3,4,5]*(1)[3,4,5,6] = ", (a2*b1).v())
print("(2)[2,3,4,5]*(3)[3,4,5,6] = ", (a2*b3).v())