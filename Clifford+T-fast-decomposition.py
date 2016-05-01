# author: Yupan Liu
# date: Apr 30, 2016
# brief: Clifford+T gate synthesis algorithm by KMM12
# based on Vadym, Maslov and Mosca's paper: arXiv:quant-ph/1206.5236v4

from sympy import *
from sympy.polys.rings import *
init_printing(use_unicode=True)

def getNorm(item):
    return re(item)**2+im(item)**2

# Calculate Initial sde
def getInitDenomExp(item):
    item = getNorm(item)
    DenomExp = 0
    DenomExpTmp = 1
    while (item*(2**DenomExpTmp)).is_Integer == False:
        print(item, DenomExp, DenomExpTmp)
        if (item*(2**(DenomExpTmp*2))).is_Integer == True:
            DenomExp += DenomExpTmp
            item *= (2**DenomExpTmp)
            DenomExpTmp = 1
        else:
            DenomExpTmp *= 2
    return DenomExp

# Get this gate's exact decomposition
def exactDecompose(U, UNow=None):
    # Hadamard gate and T gate
    Hgate = Matrix([[sqrt(Rational(1, 2)), sqrt(Rational(1, 2))],
                    [sqrt(Rational(1, 2)), -sqrt(Rational(1, 2))]], complex=True)
    Tgate = Matrix([[1, 0], [0, exp(I * pi * Rational(1, 4))]], complex=True)
    invTgate = Matrix([[1, 0], [0, exp(-I * pi * Rational(1, 4))]], complex=True)
    # Main process
    unitarySeq = ""
    sdeNow = getInitDenomExp(U[0])
    while sdeNow > 3:
        print(U, sdeNow)
        isFound = False
        invTgates = Matrix([[1,0],[0,1]], complex=True)
        invTs = ""
        for k in Range(0,4):
            if k > 0:
                invTgates = invTgates*invTgate
                invTs = invTs+"T"
            while isFound == False:
                UNow = Hgate*invTgates*U
                sdeTmp = getInitDenomExp(UNow[0])
                if sdeTmp == sdeNow-1:
                    isFound = True
                    unitarySeq = unitarySeq + invTs + "H"
                    sdeNow = sdeTmp
                    U = UNow
    unitaryRem = "Rem"
    unitarySeq = unitarySeq + unitaryRem
    return unitarySeq

# Basic quantum gate, i.e. one qubit gate of Clifford+T
Hgate = Matrix([[sqrt(Rational(1, 2)), sqrt(Rational(1, 2))],
                            [sqrt(Rational(1, 2)), -sqrt(Rational(1, 2))]], complex=True)
Sgate = Matrix([[1, 0], [0, I]], complex=True)
Tgate = Matrix([[1, 0], [0, exp(I*pi*Rational(1,4))]], complex=True)
U = Matrix([[sqrt(Rational(1,2))**10+I*sqrt(Rational(1,2))**4, (I-1)*Rational(1,2)],
            [(I+1)*Rational(1,2), sqrt(Rational(1,2))]])

print(U[0], getInitDenomExp(getNorm(U[0])), getInitDenomExp(getNorm(U[1])))
print(exactDecompose(U))