# author: Yupan Liu
# date: Nov 15, 2015
# brief: quantum circuit synthesis, from SU(d) to SU(2)
from sympy import *
init_printing(use_unicode=True)

# Input quantum gate A belongs to SU(d)
A = Matrix([[Rational(1,2), Rational(1,2), Rational(1,2), Rational(1,2)],
            [Rational(1,2), I*Rational(1,2), Rational(-1,2), I*Rational(-1,2)],
            [Rational(1,2), Rational(-1,2), Rational(1,2), Rational(-1,2)],
            [Rational(1,2), I*Rational(-1,2), Rational(-1/2), I*Rational(1/2)]])

# Calculate SU(d) operator's dimension
n = A.shape.__getitem__(0)
print(n)

# Decomposite SU(d) to two-level unitaries
U = []
for i in range(0,n-2):
    for j in range(i+1,n):
        print(i,j)
        if A[i*n+j] == 0:
            U.append(eye(n))
        else:
            # Unitary's element
            a, b = A[j], A[i*n+j]
            norm = sqrt(a*conjugate(a) + b*conjugate(b))
            # two-level unitary's SU(2) part
            Uk = zeros(n)
            Uk[i*n+i], Uk[i*n+j] = conjugate(a)/norm, conjugate(b)/norm
            Uk[j*n+i], Uk[j*n+j] = b/norm, -a/norm
            # two-level unitary's SU(2) part's inverse
            UkInv = zeros(n)
            UkInv[i*n+i], UkInv[i*n+j] = a/norm, conjugate(b)/norm
            UkInv[j*n+i], UkInv[j*n+j] = b/norm, -conjugate(a)/norm
            for k in range(0,n):
                if k != i and k != j:
                    Uk[k*n+k] = 1
            U.append(Uk)
            pprint(Uk)





