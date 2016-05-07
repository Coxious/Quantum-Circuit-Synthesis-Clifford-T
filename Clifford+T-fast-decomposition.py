# author: Yupan Liu
# date: Apr 30, 2016
# brief: Clifford+T gate synthesis algorithm by KMM12
# based on Vadym, Maslov and Mosca's paper: arXiv:quant-ph/1206.5236v4

from sympy import *
from sympy.polys.rings import *
from EuclideanDomain import EuclideanDomainNum
init_printing(use_unicode=True)

# Stack, consists of push, pop and emptyTest operations
class Stack:
    def __init__(self):
        self.items = []
    def isEmpty(self):
        return len(self.items) == 0
    def push(self, item):
        self.items.append(item)
    def pop(self):
        return self.items.pop()
    def getAll(self):
        return self.items

# Matrix representation using DenomExp enlarge
class MatrixZ:
    def __init__(self):
        rootNow = exp(pi*I*Rational(1,4))
        oneNow = AlgebraicNumber(rootNow, [0,0,0,1])
        zeroNow = AlgebraicNumber(rootNow, [0,0,0,0])
        self.m = Matrix([[oneNow, zeroNow],[zeroNow, oneNow]])
        self.de = 0

    def setValue(self, m, de):
        assert isinstance(m, Matrix)
        assert m.shape == (2,2)
        assert isinstance(m[0], AlgebraicNumber)
        assert isinstance(m[1], AlgebraicNumber)
        assert isinstance(m[2], AlgebraicNumber)
        assert isinstance(m[3], AlgebraicNumber)
        self.m = m
        self.de = de

    def __mul__(self, other):
        self.m = self.m * other.m

# Get the coefficients for (1/sqrt(2))^DenomExp w, where w in Z[1/sqrt(2),i]
# Return integer coefficients or empty (exists non-integer coefficient)
def getCoefficients(item):
    roots = [1, exp(pi*I/4), exp(pi*I/2), exp(pi*I*3/4)]
    coeff = []
    for i in range(0,len(roots)):
        terms = Stack()
        idx = 3-i
        # print("%%% Index: ", i, " RootHere: ", roots[idx])
        # Using a stack to decide whether integer coefficient exist or not
        div = expand(item / roots[idx], complex=True).as_coeff_add()
        # print("Got Div: ", div)
        for term in div:
            terms.push(term)
        # print("Initialize stack: ", terms.getAll())
        foundIntCoeff = False
        while not terms.isEmpty():
            termNow = terms.pop()
            # print("Got item: ", termNow)
            if isinstance(termNow, list) or isinstance(termNow, tuple):
                for term in termNow:
                    terms.push(term)
                # print("Update stack: ", terms.getAll())
            elif termNow.is_Integer:
                coeff.append(termNow)
                print("## Got Efficient: ", termNow)
                item -= roots[idx] * termNow
                foundIntCoeff = True
                break
        if foundIntCoeff == False:
            return []
    return coeff

# For w in Z[1/sqrt(2),i], enlarge it (i.e. multiply sqrt(2)^k) such that \tilde{w} in Z[exp(i*pi/4)]
def getDenomExp(item):
    DenomExp = 0;
    DenomExp2 = 0;
    # Enlarge it by 2
    while True:
        print("Enlarge: ", expand(item), ", DenomExp2 = ", DenomExp2)
        coeff = getCoefficients(item)
        if len(coeff) == 0:
            DenomExp2 += 1
            item *= 2
        else:
            break
    coeff = getCoefficients(item/sqrt(2))
    if len(coeff) == 0:
        DenomExp = 2*DenomExp2
    else:
        DenomExp = 2*DenomExp2-1
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
    denomExp = []
    for i in Range(0,4):
        de = getDenomExp(U[0])
        denomExp.append(de)
        U[0] *= (sqrt(2)**de)
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

# print(U[0], getInitDenomExp(getNorm(U[0])), getInitDenomExp(getNorm(U[1])))
# print(exactDecompose(U))

# testNum = (-Rational(29,2)*sqrt(2)+1222+23*I+sqrt(2)*I*Rational(37,2))*(sqrt(Rational(1,2))**89)
# testNum = sqrt(Rational(1, 2))
# de = getDenomExp(testNum)
# print(de, getCoefficients(testNum*(sqrt(2)**de)))
# print("Test Done")
