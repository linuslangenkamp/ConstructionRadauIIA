from mpmath import *
from sympy import symbols, diff, expand, Poly, solve, N, re, im

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
    roots = solve(pDer, x)

    allRoots = []
    for r in roots:
        rDps = N(r, dps)
        if rDps.is_real:
            allRoots.append(rDps)
        elif abs(im(rDps)) <= 1e-150:
            allRoots.append(re(rDps))
        else:
            rDps = N(r, 5*dps)
            allRoots.append(re(rDps))
            print("IMAG: %s" % im(rDps))

    allRoots.sort()
    allRoots = [mpf(r) for r in allRoots]

    # calculate weights
    C = [[mp.mpf(1)] * s, allRoots]

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
    outStr = f"""case IntegratorSteps::Steps{s}:
            return """
    outStr += "{"
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


### new integral approach

def fixRoots(roots, dps=150):
    allRoots = []
    for r in roots:
        rDps = N(r, dps)
        if rDps.is_real:
            allRoots.append(rDps)
        elif abs(im(rDps)) <= 1e-150:
            allRoots.append(re(rDps))
        else:
            rDps = N(r, dps)
            allRoots.append(re(rDps))

    allRoots.sort()
    return allRoots

x = symbols('x')

def lagrange(c, j, dps=150):
    p = mpf(1)
    for i in range(len(c)):
        if i != j:
            p *= (x - c[i]) / (c[j] - c[i])
    return [mpf(coeff) for coeff in Poly(expand(p)).all_coeffs()[::-1]]

def integrate(poly, roots, i, dps=150):
    S = mpf(0)
    for j in range(len(roots)):
        S += N(1, dps) / N(j+1, dps) * poly[j] * pow(roots[i], N(j+1, dps))
    return S

def generateFromIntegral(s, dps=150): 
    mp.dps = dps

    # calculate coeffs of d^(s-1) / dx^(s-1) x**(s-1) * (x-1)**s
    function = x**(s-1) * (x-1)**s
    der = expand(diff(function, x, s-1))
    pDer = Poly(der, x)
    roots = solve(pDer, x)

    roots = fixRoots(roots)

    c = [mpf(N(r, dps)) for r in roots]
    
    A = [[] for _ in range(len(c))]

    for j in range(len(c)):
        lagr = lagrange(c, j)
        for i in range(len(c)):
            a = integrate(lagr, c, i)
            A[i].append(a)
    return A, inverse(A), c

import time

for m in range(2, 50):
    startTime = time.time()
    A, Ainv, c = generateFromIntegral(m)
    print(f"{m} {time.time() - startTime}")
"""
for s in range(1, 37):
    genCppCode(s, 150, 53)
    print("")
"""

