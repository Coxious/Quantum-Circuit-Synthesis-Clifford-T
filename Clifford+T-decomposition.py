# author: Yupan Liu
# date: Nov 16, 2015
# brief: Clifford+T gate synthesis algorithm by Matsumoto and Amano's normal form
# based on Giles and Selinger's paper: arXiv quant-ph/1312.6584v1
from sympy import *
from sympy.polys.rings import *
init_printing(use_unicode=True)

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
    return expand(UBloch, complex=True)

# Calculate the Bloch rep's least denominator exponent
def getLeastDenominatorExponent(U):
    return []

# Basic quantum gate, i.e. one qubit gate of Clifford+T
Hgate = Matrix([[sqrt(Rational(1,2)), sqrt(Rational(1,2))],
                [sqrt(Rational(1,2)), -sqrt(Rational(1,2))]], complex=True)
Sgate = Matrix([[1, 0], [0, I]], complex=True)
Tgate = Matrix([[1, 0], [0, exp(I*pi*Rational(1,4))]], complex=True)
U = Matrix([[I,sqrt(2)],[-sqrt(2),I]])

pprint(Hgate)
pprint(Sgate)
pprint(Tgate)
pprint(getBlochRep(Hgate))
pprint(getBlochRep(Sgate))
pprint(getBlochRep(Tgate))
pprint(getBlochRep(U))
