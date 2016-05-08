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
    def __init__(self, m):
        self.root = exp(pi*I*Rational(1,4))
        assert isinstance(m, list)
        assert len(m) == 4
        assert isinstance(m[0], EuclideanDomainNum)
        assert isinstance(m[1], EuclideanDomainNum)
        assert isinstance(m[2], EuclideanDomainNum)
        assert isinstance(m[3], EuclideanDomainNum)
        self.m = m

    def __getitem__(self, index):
        assert isinstance(index, int) or isinstance(index, Integer)
        assert 0 <= index <= 3
        return self.m[index]

    def __setitem__(self, index, value):
        assert isinstance(index, int) or isinstance(index, Integer)
        assert isinstance(value, EuclideanDomainNum)
        assert 0 <= index <= 3
        self.m[index] = value

    def __add__(self, other):
        m = [0, 0, 0, 0]
        for i in range(0, 4):
            m[i] = self.m[i] + other.m[i]
        return MatrixZ(m)

    def __sub__(self, other):
        m = [0, 0, 0, 0]
        for i in range(0, 4):
            m[i] = self.m[i] - other.m[i]
        return MatrixZ(m)

    def __mul__(self, other):
        # print("Multiplication Running!")
        r, s = self.m, other.m
        # print(r[0].v(), s[0].v(), (r[0]*s[0]).v(), (r[0]*s[0]).getDenomExp())
        # print(r[1].v(), s[2].v(), (r[1]*s[2]).v(), (r[1]*s[2]).getDenomExp())
        # print((r[0]*s[0]+r[1]*s[2]).v())
        m = [r[0]*s[0]+r[1]*s[2],\
             r[0]*s[1]+r[1]*s[3],\
             r[2]*s[0]+r[3]*s[2],\
             r[2]*s[1]+r[3]*s[3]]
        # print("Multiplication Done")
        return MatrixZ(m)

    def getMatrix(self):
        return Matrix([[self.m[0].v(), self.m[1].v()], [self.m[2].v(), self.m[3].v()]], complex=True)

# For w in Z[1/sqrt(2),i], enlarge it (i.e. multiply sqrt(2)^k) such that \tilde{w} in Z[exp(i*pi/4)]
def getDenomExp(item):
    DenomExp = 0;
    DenomExp2 = 0;
    # Enlarge it by 2
    while True:
        # print("Enlarge: ", expand(item), ", DenomExp2 = ", DenomExp2)
        coeff = getCoefficients(item)
        if len(coeff) == 0:
            DenomExp2 += 1
            item *= 2
        else:
            break
    coeff = getCoefficients(item / sqrt(2))
    if len(coeff) == 0:
        DenomExp = 2 * DenomExp2
    else:
        DenomExp = 2 * DenomExp2 - 1
    return DenomExp

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
                # print("## Got Efficient: ", termNow)
                item -= roots[idx] * termNow
                foundIntCoeff = True
                break
        if foundIntCoeff == False:
            return []
    return coeff

# Get this gate's exact decomposition
def exactDecompose(UOri, UNow=None):
    # Hadamard gate and T gate
    one = EuclideanDomainNum([0,0,0,1],0)
    zero = EuclideanDomainNum([0,0,0,0],0)
    Hgate = MatrixZ([EuclideanDomainNum([0,0,0,1],1), \
                    EuclideanDomainNum([0,0,0,1],1),\
                    EuclideanDomainNum([0,0,0,1],1), \
                    EuclideanDomainNum([0,0,0,-1],1)])
    Tgate = MatrixZ([one, zero, zero, EuclideanDomainNum([0,0,1,0],0)])
    invTgate = MatrixZ([one, zero, zero, EuclideanDomainNum([-1,0,0,0],0)])
    U = MatrixZ([one, one, one, one])
    # Main process
    unitarySeq = ""
    denomExp = []
    for i in Range(0,4):
        de = getDenomExp(UOri[i])
        denomExp.append(de)
        U[i] = EuclideanDomainNum(getCoefficients(UOri[i]*(sqrt(2)**de)),de)
    sdeNow = U[0].getDenomExp()*4 - U[0].gdeForNorm()
    while sdeNow > 3:
        print(U.getMatrix(), sdeNow)
        isFound = False
        invTgates = MatrixZ([one, zero, zero, one])
        invTs = ""
        for k in Range(0,4):
            print("Running: ", k)
            if k > 0:
                invTgates = invTgates*invTgate
                invTs = invTs+"T"
            print("Running Mul Done: ", k)
            UNow = U
            while isFound == False:
                print("Running While: ", k)
                UNow = Hgate*invTgates*UNow
                print("Running UNow Done: ", k)
                sdeTmp = UNow[0].getDenomExp()*4 - UNow[0].gdeForNorm()
                print("Running sdeTmp Done: ", k)
                print(UNow.getMatrix(), sdeTmp)
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

print(exactDecompose(U))

testCoeffAnalysis = False
if testCoeffAnalysis:
    testNum = (-Rational(29, 2) * sqrt(2) + 1222 + 23 * I + sqrt(2) * I * Rational(37, 2)) * (
    sqrt(Rational(1, 2)) ** 89)
    testNum1 = sqrt(Rational(1, 2))
    testNum2 = -sqrt(Rational(1, 2))
    testNum3 = 1
    testNum4 = 0
    testNum5 = exp(I * pi * Rational(1, 4))
    testNum6 = exp(-I * pi * Rational(1, 4))
    # de = getDenomExp(testNum)
    # print("Test0: ", de, getCoefficients(testNum*(sqrt(2)**de)))
    de = getDenomExp(testNum1)
    print("Test1: ", de, getCoefficients(testNum1*(sqrt(2)**de)))
    de = getDenomExp(testNum2)
    print("Test2: ", de, getCoefficients(testNum2*(sqrt(2)**de)))
    de = getDenomExp(testNum3)
    print("Test3: ", de, getCoefficients(testNum3*(sqrt(2)**de)))
    de = getDenomExp(testNum4)
    print("Test4: ", de, getCoefficients(testNum4*(sqrt(2)**de)))
    de = getDenomExp(testNum5)
    print("Test5: ", de, getCoefficients(testNum5*(sqrt(2)**de)))
    de = getDenomExp(testNum6)
    print("Test6: ", de, getCoefficients(testNum6*(sqrt(2)**de)))
print("Test Done")
