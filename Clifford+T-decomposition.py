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
    return expand(UBloch, complex=True, rational=True)


# Calculate the Bloch rep's least denominator exponent
def getLeastDenominatorExponent(U):
    leastDenoExp = 0
    for u in U:
        while True:
            # At most k times
            val = expand(u*(sqrt(2)**leastDenoExp), complex=True)
            x, y = val.as_coeff_Add()
            y = y/sqrt(2)
            if x.is_Integer and y.is_Integer:
                break
            else:
                # Maintain the global minimum of its least denominator exponent
                leastDenoExp = leastDenoExp+1
    return leastDenoExp


# Calculate the Bloch rep's k-parity
def getKParity(U, k):
    #pprint(expand(U*(sqrt(2)**k), complex=True))
    UParity = zeros(3)
    UParityLoc = 0
    for u in U:
        val = expand(u*(sqrt(2)**k), complex=True)
        x, y = val.as_coeff_Add()
        UParity[UParityLoc] = x%2
        UParityLoc += 1
    return UParity


# Get its Matsumoto-Amano normal form
def getMAFormType(U):
    #pprint(U)
    if U==Matrix([[1,1,0], [1,1,0], [0,0,0]]):
        return 1;
    elif U==Matrix([[0,0,0], [1,1,0], [1,1,0]]):
        return 2;
    elif U==Matrix([[1,1,0], [0,0,0], [1,1,0]]):
        return 3;


# Get this gate's exact decomposition
def exactDecompose(U):
    unitarySeq = ""
    LeastDenoExp = getLeastDenominatorExponent(U)
    while(LeastDenoExp >= 0):
        UType = 0
        UParity = getKParity(U, LeastDenoExp)
        pprint(UParity)
        if LeastDenoExp > 0:
            UType = getMAFormType(UParity)
        if UType == 0 and UParity != DiagonalMatrix(1,1,1):
            unitarySeq += "C"
        elif UType == 1:
            unitarySeq += "T"
            U = TBloch*U
        elif UType == 2:
            unitarySeq += "HT"
            U = TBloch*HBloch*U
        elif UType == 3:
            unitarySeq += "SHT"
            U = TBloch*HBloch*SBloch*U
        LeastDenoExp -= 1
    return unitarySeq


# Basic quantum gate, i.e. one qubit gate of Clifford+T
Hgate = Matrix([[sqrt(Rational(1,2)), sqrt(Rational(1,2))],
                [sqrt(Rational(1,2)), -sqrt(Rational(1,2))]], complex=True)
Sgate = Matrix([[1, 0], [0, I]], complex=True)
Tgate = Matrix([[1, 0], [0, exp(I*pi*Rational(1,4))]], complex=True)
U = Matrix([[sqrt(Rational(1,2)), (I-1)*Rational(1,2)],
            [(I+1)*Rational(1,2), sqrt(Rational(1,2))]])

# Clifford+T gate
HBloch = getBlochRep(Hgate)
SBloch = getBlochRep(Sgate)
TBloch = getBlochRep(Tgate)
UBloch = getBlochRep(Tgate*Hgate*Tgate*Sgate*Hgate*Tgate*Sgate*Hgate*Tgate*Hgate*Tgate*Hgate*Tgate*Sgate*Sgate)
pprint(UBloch)
pprint(getLeastDenominatorExponent(UBloch))
print(exactDecompose(UBloch))

# Leave problem:
# 1. Clifford Gate at the tail -> Decomposition?
# 2. Clifford Gate -> Permutation -> Affect on MA Form?
# 3. Tooooo Slow!