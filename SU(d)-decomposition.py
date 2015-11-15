# author: Yupan Liu
# date: Nov 15, 2015
# brief: quantum circuit synthesis, from SU(d) to SU(2)
from sympy import *
init_printing(use_unicode=True)

# Get the two-level unitary's non-trivial basis
def getBaseKet(ket, qubitCnt):
    baseKet = ""
    for i in range(qubitCnt):
        if ket&(2**i):
            baseKet = "1" + baseKet
        else:
            baseKet = "0" + baseKet
    return baseKet

# Generate Grey Code from ket |s> to ket |t>
def generateGreyCode(baseKetS, baseKetT, n):
    GreyBasis = []
    GreyBasis.append([baseKetS,-1])
    # Get every Grey code between s and t
    for i in range(n):
        if baseKetS[i] != baseKetT[i]:
            if baseKetS[i] == "0":
                baseKetS = baseKetS[0:i]+"1"+baseKetS[i+1:n]
            else:
                baseKetS = baseKetS[0:i]+"0"+baseKetS[i+1:n]
            GreyBasis.append([baseKetS,i])
    return GreyBasis

# Generate SU(2) Gates by CNOT and single qubit gate
def generateCNOTs(baseKet, targetLoc, n):
    CNOTseq = []
    return CNOTseq

# Approximate SU(2) quantum gates
def SU2Approx(twoLevelUnitary, s, t, n):
    U = zeros(2)
    U[0] = twoLevelUnitary[s*n+s]
    U[1] = twoLevelUnitary[s*n+t]
    U[2] = twoLevelUnitary[t*n+s]
    U[3] = twoLevelUnitary[t*n+t]
    return U

# Input quantum gate A belongs to SU(d)
# Examples: 2 qubit quantum Fourier transform
A = Matrix([[Rational(1,2), Rational(1,2), Rational(1,2), Rational(1,2)],
            [Rational(1,2), I*Rational(1,2), Rational(-1,2), I*Rational(-1,2)],
            [Rational(1,2), Rational(-1,2), Rational(1,2), Rational(-1,2)],
            [Rational(1,2), I*Rational(-1,2), Rational(-1/2), I*Rational(1/2)]])

# Calculate SU(d) operator's dimension
n = A.shape.__getitem__(0)
print(n)

# Decompose SU(d) to two-level unitaries
U = []
Uidx = []
for j in range(0,n-1):
    for i in range(j+1,n):
        # generate the i*n+j-th two-level unitary
        # print(i,j)
        Uidx.append([j,i])
        if A[i*n+j] == 0:
            U.append(eye(n))
        else:
            # The last two-level unitary maybe is not in SU(2)
            if j == n-2 and i == n-1:
                U.append(expand(A))
                continue
            # Unitary's element
            a, b = A[j*n+j], A[i*n+j]
            norm = sqrt(a*conjugate(a) + b*conjugate(b))
            # two-level unitary's SU(2) part
            UkInv = zeros(n)
            UkInv[j*n+j], UkInv[j*n+i] = conjugate(a)/norm, conjugate(b)/norm
            UkInv[i*n+j], UkInv[i*n+i] = b/norm, -a/norm
            # two-level unitary's SU(2) part's inverse
            Uk = zeros(n)
            Uk[j*n+j], Uk[j*n+i] = a/norm, conjugate(b)/norm
            Uk[i*n+j], Uk[i*n+i] = b/norm, -conjugate(a)/norm
            for k in range(0,n):
                if k != i and k != j:
                    Uk[k*n+k] = UkInv[k*n+k] = 1
            # The k-th two-level unitary
            U.append(expand(Uk))
            # pprint(expand(UkInv))
            # Update the matrix A to guarantee A[i,j] = 0
            A = UkInv*A
            # pprint(expand(A))
# verify whether the decomposition is correct or not.
V = eye(n)
for Ui in U:
    pprint(Ui)
    V = V*Ui
pprint(expand(V))

# Transform two-level unitary to CNOT and single qubit gates
idx = 0
SU2Gates = []
for Ui in U:
    s, t = Uidx[idx][0], Uidx[idx][1]
    idx = idx+1
    # Generate Grey Codes
    BaseS, BaseT = getBaseKet(s,n), getBaseKet(t,n)
    GreyBasisSeq = generateGreyCode(BaseS, BaseT, n)
    print(GreyBasisSeq)
    GreyBasisLength = len(GreyBasisSeq)
    # Decompose the two-level unitary, i.e. permute the unitary to get single qubit gate
    for i in range(1,GreyBasisLength):
        SU2Gates.append(generateCNOTs(GreyBasisSeq[i][0], GreyBasisSeq[i][1], n))
    # Approximate single qubit gate
    SU2Gates.append(SU2Approx(Ui, s, t, n))
    # Decompose the two-level unitary, i.e. permute the unitary to restore original version
    for i in range(GreyBasisLength-1,-1,-1):
        if i < 1:
            break
        SU2Gates.append(generateCNOTs(GreyBasisSeq[i][0], GreyBasisSeq[i][1], n))
