# author: Yupan Liu
# date: Nov 16, 2015
# brief: single qubit gate approximation by Z[i,1/sqrt(2)] unitary
# based on Peter Selinger's paper: arXiv quant-ph/1212.6253v2
from sympy import *

epsilon = Rational(1/50)
C = Rational(5,2)+2*log(1+sqrt(2))/log(2)
# k is Unitary's least denominator exponent
k = ceiling(C+2*log(1/epsilon)/log(2))
#
n = floor(4*sqrt(2)/epsilon)


# generate new u\Bar randomly
def generateUBar():
    return 0


# solve Diophantine Equation t^{\dagger}t = \xi
def solveDiophantineEq(xi):
    return 0


# z-rotation operator approximation
def zRotation(theta):
    U = zeros(2)
    while 1:
        uBar = generateUBar()
        u = uBar*(sqrt(2)**k)
        xi = 2**k-conjugate(u)*u
        t = solveDiophantineEq(xi)
        if t == -1:
            uBar = generateUBar()
        else:
            U[0], U[1], U[2], U[3] = u, -conjugate(t), t, conjugate(t)
            U = expand((sqrt(Rational(1,2))**k)*U)
            break
    pprint(U)
    return U


# get unitary's Enler angle
def getEulerAngle(U):
    return []

a = I
b = 1
U = [[a, b], [conjugate(b), -conjugate(a)]]
[alpha, beta, gamma] = getEulerAngle(U)
pprint(zRotation(alpha))
pprint(zRotation(beta))
pprint(zRotation(gamma))