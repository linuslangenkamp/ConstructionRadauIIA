from mpmath import *
from sympy import symbols, diff, expand, Poly, solve, N, re, im


def genCppCode(s, dps=150, mindps=53):
    A, A_inv, c, invRowSum = generateFromIntegral(s, dps)
    mp.dps = mindps
    outStr = f"""case IntegratorSteps::Steps{s}:
            return """
    outStr += "{"

    # c vector
    outStr += "{"
    for i in range(s):
        outStr += str(c[i]) + ",\n"
    outStr = outStr[:-2] +"},\n"

    # b vector
    outStr += "{"
    for i in range(s):
        outStr += str(A[-1][i]) + ",\n"
    outStr = outStr[:-2] +"},\n"

    # ainv matrix
    outStr += "{"
    for i in range(s):
        outStr += "{"
        for j in range(s):
            outStr += str(A_inv[i, j]) + ",\n"
        outStr = outStr[:-2] + "},\n"
    outStr = outStr[:-2] + "},\n"

    # row sum of ainv
    outStr += "{"
    for i in range(s):
        outStr += str(invRowSum[i]) + ",\n"
    outStr = outStr[:-2] + "},\n"

    # step count
    outStr += str(s)

    outStr += "};"
    return outStr


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
    inv = inverse(A)
    rowSum = [sum(row) for row in inv.tolist()]
    return A, inverse(A), c, rowSum

import time

for m in range(2, 37):
    with open("radau1_37.txt", "w") as f:
        startTime = time.time()
        #f.write()
        print(genCppCode(m))
        #print(f"{m} {time.time() - startTime}")
