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

# Get the Bloch sphere representation
def getBlochRep(U):
    # Renormalize the operator from U(2) to SU(2)
    phaseFactor = sqrt(det(U))
    a, b = expand(U[0]*phaseFactor), expand(U[1]*phaseFactor)
    # Calculate its Bloch sphere representation
    UBloch = zeros(3)
    UBloch[0] = Rational(1,2)*(pow(a,2)+pow(conjugate(a),2)-pow(b,2)-pow(conjugate(b),2))
    UBloch[1] = I*Rational(1,2)*(pow(conjugate(a),2)+pow(conjugate(b),2)-pow(a,2)-pow(b,2))
    UBloch[2] = -conjugate(a)*conjugate(b)-a*b
    UBloch[3] = I*Rational(1,2)*(pow(a,2)+pow(conjugate(b),2)-pow(conjugate(a),2)-pow(b,2))
    UBloch[4] = Rational(1,2)*(pow(conjugate(a),2)+pow(conjugate(b),2)+pow(a,2)+pow(b,2))
    UBloch[5] = I*(conjugate(a)*conjugate(b)-a*b)
    UBloch[6] = conjugate(a)*b+conjugate(b)*a
    UBloch[7] = I*(conjugate(a)*b-conjugate(b)*a)
    UBloch[8] = conjugate(a)*a-conjugate(b)*b
    UBloch.transpose()
    # Guarantee the result is simple enough
    return expand(UBloch, complex=True, rational=True)

# generate new u\Bar randomly
def generateUBar():
    return 0


# solve Diophantine Equation t^{\dagger}t = \xi
def solveDiophantineEq(xi):
    # xi = x + y sqrt(2)
    x, y = xi.as_coeff_Add()
    y = y/sqrt(2)
    # find h such that p|h^2+1
    p = x**2 - 2*y**2
    b = 1
    for i in range(1, p):
        if i**(Rational(p-1,2))%p == p-1:
            b = i
            break
    h = b**Rational(p-1,4)
    # gcd() on Z[i,sqrt(i)]
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
    R = getBlochRep(U)
    alpha = atan2(R[5], R[8])
    beta = atan2(-R[2], sqrt(R[5]**2+R[8]**2))
    gamma = atan2(R[1], R[0])
    if alpha == nan:
        alpha = 0
    if beta == nan:
        beta = 0
    if gamma == nan:
        gamma = 0
    return [alpha, beta, gamma]

a = I
b = 1
U = [[a, b], [conjugate(b), -conjugate(a)]]
[alpha, beta, gamma] = getEulerAngle(U)
# only z-rotation can be approximated
pprint(zRotation(alpha))
pprint(zRotation(beta))
pprint(zRotation(gamma))