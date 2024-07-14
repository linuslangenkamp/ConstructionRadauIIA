from mpmath import *
from sympy import symbols, diff, expand, Poly

def generateRadauIIA(s, dps, mindps, printout=False):
    """
    Generates Radau IIA s step schema with dps digits and returns value up to mindps digits
    """
    mp.dps = dps
    x = symbols('x')

    # calculate coeffs of d^(s-1) / dx^(s-1) x**(s-1) * (x-1)**s
    function = x**(s-1) * (x-1)**s
    der = expand(diff(function, x, s-1))
    pDer = Poly(der, x)
    coefficients = [int(pDer.coeff_monomial(x**(s-k))) for k in range(s+1)]

    # calculate weights
    C = [[mp.mpf(1)] * s, polyroots(coefficients)]
    rhs = [C[1]]
    for k in range(2, s+1):
        C.append([C[1][j]**k for j in range(s)])
        rhs.append([C[k][j] / k for j in range(s)])

    # construct matrix
    zero = mp.mpf(0)
    matrixC = []
    for i in range(s):
        zeroR = [zero] * s
        for j in range(s):
            matrixC.append(zeroR * j + C[i] + zeroR * (s-1-j))

    # flatten rhs
    rhsFlat = [rhs[i][j] for i in range(s) for j in range(s)]

    # solve
    x = lu_solve(matrixC, rhsFlat)

    # construct a_ij
    A = []
    for i in range(s):
        row = []
        for j in range(s):
            row.append(x[s * i + j])
        A.append(row)

    # calculate inverse of A
    A_inv = inverse(A)

    invRowSum = []
    for i in range(s):
        S = sum(A_inv[i, j] for j in range(s))
        invRowSum.append(S)

    mp.dps = mindps
    if printout:
        # print scheme
        print("KNODS:")
        for i in range(s):
            print("c[%s] = %s" % (i+1, C[1][i]))

        print("COEFFS:")
        for i in range(s):
            for j in range(s):
                print("A[%s, %s] = %s" % (i + 1, j + 1, x[s * i + j]))

        print("INVERSE COEFFS:")
        for i in range(s):
            for j in range(s):
                print("inv(A)[%s, %s] = %s" % (i + 1, j + 1, A_inv[i, j]))
    return A, A_inv, C[1], invRowSum


def genCppCode(s, dps, mindps):
    A, A_inv, c, invRowSum = generateRadauIIA(s, dps, mindps)
    mp.dps = mindps
    outStr = "{"
    outStr += "{"
    for i in range(s):
        outStr += str(c[i]) + ",\n"
    outStr = outStr[:-2] +"},\n"

    outStr += "{"
    for i in range(s):
        outStr += "{"
        for j in range(s):
            outStr += str(A[i][j]) + ",\n"
        outStr = outStr[:-2] +"},\n"
    outStr = outStr[:-2] + "},\n"

    outStr += "{"
    for i in range(s):
        outStr += "{"
        for j in range(s):
            outStr += str(A_inv[i, j]) + ",\n"
        outStr = outStr[:-2] + "},\n"
    outStr = outStr[:-2] + "},\n"

    outStr += "{"
    for i in range(s):
        outStr += str(invRowSum[i]) + ",\n"
    outStr = outStr[:-2] + "},\n"

    outStr += str(s)

    outStr += "};"
    print(outStr)
    return True

for s in range(1, 8):
    print("Schema: %s" % s)
    genCppCode(s, 150, 53)
    print("")