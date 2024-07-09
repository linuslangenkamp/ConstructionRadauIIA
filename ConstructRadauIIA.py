from mpmath import *
from sympy import symbols, diff, expand, Poly

s = 4
x = symbols('x')

function = x**(s-1) * (x-1)**s
der = expand(diff(function, x, s-1))
pDer = Poly(der, x)
coefficients = [int(pDer.coeff_monomial(x**(s-k))) for k in range(s+1)]
mp.dps = 100

C = [[mp.mpf(1)] * s, polyroots(coefficients)]
rhs = [C[1]]
for k in range(2, s+1):
    C.append([C[1][j]**k for j in range(s)])
    rhs.append([C[k][j] / k for j in range(s)])

zero = mp.mpf(0)
Cmatrix = []

for i in range(s):
    rowR = C[i]
    zeroR = [zero] * s
    for j in range(s):
        Cmatrix.append(zeroR * j + C[i] + zeroR * (s-1-j))

rhsflattended = [rhs[i][j] for i in range(s) for j in range(s)]

x = lu_solve(Cmatrix, rhsflattended)

A = []
for i in range(s):
    row = []
    for j in range(s):
        row.append(x[s * i + j])
    A.append(row)
A_inv = inverse(A)

print("KNODS:")
for i in range(s):
    print("c[%s] = %s" % (i+1, C[1][i]))

print("COEFFS:")
for i in range(s):
    row = []
    for j in range(s):
        print("A[%s, %s] = %s" % (i + 1, j + 1, x[s * i + j]))

print("INVERSE COEFFS:")

for i in range(s):
    for j in range(s):
        print("inv(A)[%s, %s] = %s" % (i + 1, j + 1, A_inv[i, j]))