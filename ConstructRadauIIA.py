from mpmath import *
from sympy import symbols, diff, expand, Poly, solve, N, re, im
import time


def genCppCode(s, dps=150, mindps=40):
    c, c0, cBisection, c0Bisection, b, D1, D2, LB, LB0 = generate(s, dps)
    mp.dps = mindps

    outStr = f"""case IntegratorSteps::Steps{s}:
            return """
    outStr += "{"

    # c
    outStr += "{"
    for i in range(len(c)):
        outStr += str(c[i]) + ","
    outStr = outStr[:-1] +"},\n"

    # c0
    outStr += "{"
    for i in range(len(c0)):
        outStr += str(c0[i]) + ","
    outStr = outStr[:-1] +"},\n"

    # cBisection
    outStr += "{"
    for i in range(len(cBisection)):
        outStr += str(cBisection[i]) + ","
    outStr = outStr[:-1] +"},\n"

    # c0Bisection
    outStr += "{"
    for i in range(len(c0Bisection)):
        outStr += str(c0Bisection[i]) + ","
    outStr = outStr[:-1] +"},\n"

    # b
    outStr += "{"
    for i in range(len(b)):
        outStr += str(b[i]) + ","
    outStr = outStr[:-1] +"},\n"

    # D1 matrix
    outStr += "{"
    for i in range(len(D1)):
        outStr += "{"
        for j in range(len(D1[0])):
            outStr += str(D1[i][j]) + ","
        outStr = outStr[:-1] + "},\n"
    outStr = outStr[:-2] + "},\n"

    # D2 matrix
    outStr += "{"
    for i in range(len(D2)):
        outStr += "{"
        for j in range(len(D2[0])):
            outStr += str(D2[i][j]) + ","
        outStr = outStr[:-2] + "},\n"
    outStr = outStr[:-2] + "},\n"

    # LB matrix
    outStr += "{"
    for i in range(len(LB)):
        outStr += "{"
        for j in range(len(LB[0])):
            outStr += str(LB[i][j]) + ","
        outStr = outStr[:-2] + "},\n"
    outStr = outStr[:-2] + "},\n"

    # LB0 matrix
    outStr += "{"
    for i in range(len(LB0)):
        outStr += "{"
        for j in range(len(LB0[0])):
            outStr += str(LB0[i][j]) + ","
        outStr = outStr[:-1] + "},\n"
    outStr = outStr[:-2] + "},\n"

    # step count
    outStr += str(s) + "};\n"

    return outStr


x = symbols('x')

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

def generate(s, dps=150): 
    mp.dps = dps

    # calculate coeffs of d^(s-1) / dx^(s-1) x**(s-1) * (x-1)**s
    function = x**(s-1) * (x-1)**s
    der = expand(diff(function, x, s-1))
    pDer = Poly(der, x)
    roots = solve(pDer, x)
    roots = fixRoots(roots)

    # radau nodes
    c = [mpf(N(r, dps)) for r in roots]

    # radau nodes with 0
    c0 =  [mpf(0)] + c

    # radau nodes after bisection
    cBisection = [mpf(1/2 * k) + mpf(N(r/2, dps)) for k in range(2) for r in roots]

    # radau nodes after bisection with 0
    c0Bisection = [mpf(0)] + [mpf(1/2 * k) + mpf(N(r/2, dps)) for k in range(2) for r in roots]

    # quadrature weights
    b = []

    # helper for derivative matrices D1, D2
    weights = [weight(c0, j) for j in range(len(c0))]
    D1 = [[None for _ in range(len(c0))] for _ in range(len(c0))] # diff matrix at c0
    D2 = [[None for _ in range(len(c0))] for _ in range(len(c0))] # diff 2 matrix at c0

    # lagrange evaluation at bisection points
    LB = []
    LB0 = []

    # eval d/dx L
    for i in range(s + 1):
        lagr = lagrange(c0, i)
        if i > 0:
            if s == 1:
                b = [mpf(1)]
            else:
                b.append(integrate(lagr, c0, s))
        for j in range(s + 1):
            if i != j:
                D1[i][j] = weights[j] / (weights[i] * (c0[i] - c0[j]))
        D1[i][i] = -sum(D1[i][j] for j in range(s + 1) if j != i)
    
    # eval d^2/dx^2 L
    for i in range(s + 1):
        for j in range(s + 1):
            if i != j:
                D2[i][j] = D1[i][j] * (D1[i][i] - D1[j][j]) / (c0[i] - c0[j])
        D2[i][i] = -sum(D2[i][j] for j in range(s + 1) if j != i)

    # lagrange evaluation at cBisection
    for coll in cBisection:
        lagrC = []
        for k in range(s):
            factor = mpf(1)
            for d in range(s):
                if (d != k):
                    factor *= (coll - c[d]) / (c[k] - c[d])
            lagrC.append(factor)
        LB.append(lagrC)

    for coll in c0Bisection:
        lagrC = []
        for k in range(s + 1):
            factor = mpf(1)
            for d in range(s + 1):
                if (d != k):
                    factor *= (coll - c0[d]) / (c0[k] - c0[d])
            lagrC.append(factor)
        LB0.append(lagrC)
    return c, c0, cBisection, c0Bisection, b, D1, D2, LB, LB0


for m in range(1, 71):
    with open("radauNew.txt", "a") as f:
        startTime = time.time()
        f.write(genCppCode(m))
        print(f"{m} {time.time() - startTime}")
