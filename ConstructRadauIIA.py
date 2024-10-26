from mpmath import *
from sympy import symbols, diff, expand, Poly, solve, N, re, im


def genCppCode(s, dps=150, mindps=53):
    A, A_inv, c, invRowSum = generateFromDerivatives(s, dps)
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
    return p

def lagrangeToCoeffs(p, dps=150):
    return [mpf(coeff) for coeff in Poly(expand(p)).all_coeffs()[::-1]]

def integrate(poly, roots, i, dps=150):
    poly = lagrangeToCoeffs(poly, dps=dps)
    S = mpf(0)
    for j in range(len(roots)):
        S += N(1, dps) / N(j+1, dps) * poly[j] * pow(roots[i], N(j+1, dps))
    return S

def weight(c, i, dps=150):
    mp.dps = dps
    p = mpf(1)
    for j in range(len(c)):
        if i != j:
            p *= (c[i] - c[j])
    return mpf(1) / p

def generateFromDerivatives(s, dps=150): 
    mp.dps = dps

    # calculate coeffs of d^(s-1) / dx^(s-1) x**(s-1) * (x-1)**s
    function = x**(s-1) * (x-1)**s
    der = expand(diff(function, x, s-1))
    pDer = Poly(der, x)
    roots = solve(pDer, x)

    roots = fixRoots(roots)

    c = [mpf(N(r, dps)) for r in roots]
    c0 =  [mpf(0)] + c
    cBisection = [mpf(1/2 * k) + mpf(N(r/2, dps)) for k in range(2) for r in roots]
    c0Bisection = [mpf(0)] + [mpf(1/2 * k) + mpf(N(r/2, dps)) for k in range(2) for r in roots]
    b = []
    weights = [weight(c0, j) for j in range(len(c0))]
    D1 = [[None for _ in range(len(c0))] for _ in range(len(c0))] # diff matrix at c0
    D2 = [[None for _ in range(len(c0))] for _ in range(len(c0))] # diff 2 matrix at c0

    LB = [[None for _ in range(len(c0))] for _ in range(len(c0))] # lagrange poly basis evaluated at cBisection
    # eigentlich brauche ich nur L mit c als basis knoten -> eval at 0 -> normale interpolation
    # überlege was ich hier benötige
    for i in range(s + 1):
        lagr = lagrange(c0, i)
        if i > 0:
            b.append(integrate(lagr, c0, s))
        for j in range(s + 1):
            if i != j:
                D1[i][j] = weights[j] / (weights[i] * (c0[i] - c0[j]))
        D1[i][i] = -sum(D1[i][j] for j in range(s + 1) if j != i)

    for i in range(s + 1):
        for j in range(s + 1):
            if i != j:
                D2[i][j] = D1[i][j] * (D1[i][i] - D1[j][j]) / (c0[i] - c0[j])
        D2[i][i] = -sum(D2[i][j] for j in range(s + 1) if j != i)
    return


import time

for m in range(3, 61):
    with open("radau1_37.txt", "w") as f:
        startTime = time.time()
        #f.write()
        print(generateFromDerivatives(m))
        print(f"{m} {time.time() - startTime}")
