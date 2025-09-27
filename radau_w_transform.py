from mpmath import *
from sympy import symbols, diff, expand, Poly, solve, N, re, im
import time

def genCppCode(clist, c0list, blist, Dlist, w0list, wlist, dps=250, mindps=40):
    mp.dps = mindps

    with open("radauConstantsC.txt", "w") as f:
        # c
        outStr = "{"
        for c in clist:
            for i in range(len(c)):
                outStr += str(c[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # c0
    with open("radauConstantsC0.txt", "w") as f:
        outStr = "{"
        for c0 in c0list:
            for i in range(len(c0)):
                outStr += str(c0[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # b
    with open("radauConstantsB.txt", "w") as f:
        outStr = "{"
        for b in blist:
            for i in range(len(b)):
                outStr += str(b[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # D1 matrix
    with open("radauConstantsD.txt", "w") as f:
        outStr = "{"
        for D1 in Dlist:
            for i in range(len(D1)):
                for j in range(len(D1[0])):
                    outStr += str(D1[i][j]) + ","
                outStr += "\n"
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # barycentric weights (incl. 0)
    with open("radauConstantsW0.txt", "w") as f:
        outStr = "{"
        for w0 in w0list:
            for i in range(len(w0)):
                outStr += str(w0[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    # barycentric weights (excl. 0)
    with open("radauConstantsW.txt", "w") as f:
        outStr = "{"
        for w in wlist:
            for i in range(len(w)):
                outStr += str(w[i]) + ","
            outStr += "\n"
        outStr = outStr[:-2] + "}\n"
        f.write(outStr)

    return outStr


x = symbols("x")


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
        S += (
            N(1, dps)
            / N(j + 1, dps)
            * poly[j]
            * (pow(N(1, dps), j + 1) - pow(N(0, dps), j + 1))
        )
    return S


def weight(c, i, dps=150):
    mp.dps = dps
    p = mpf(1)
    for j in range(len(c)):
        if i != j:
            p *= c[i] - c[j]
    return mpf(1) / p


def tridiagonal_eigenvalues(n, dps=150):
    mp.dps = dps

    alpha = mp.mpf(1)

    # tridiag entries
    a_0 = -mp.mpf(1) / mp.mpf(3)
    a_j = lambda j: -mp.mpf(1) / (
        (mp.mpf(2) * j + mp.mpf(1)) * (mp.mpf(2) * j + mp.mpf(3))
    )

    b_1 = mp.sqrt(mp.mpf(8) / (mp.mpf(9) * (mp.mpf(3) + alpha)))
    b_j = lambda j: mp.sqrt(
        (mp.mpf(4) * j * j * (j + mp.mpf(1)) * (j + mp.mpf(1)))
        / (
            (mp.mpf(2) * j)
            * (mp.mpf(2) * j + mp.mpf(1))
            * (mp.mpf(2) * j + mp.mpf(1))
            * (mp.mpf(2) * j + mp.mpf(2))
        )
    )

    # set matrix
    matrix = mp.zeros(n)
    for i in range(n):
        if i == 0:
            matrix[i, i] = a_0  # First diagonal element
            if n > 1:
                matrix[i, i + 1] = b_1
                matrix[i + 1, i] = b_1
        else:
            matrix[i, i] = a_j(i)  # Diagonal element
            if i < n - 1:  # Off-diagonal elements
                b_val = b_j(i + 1)
                matrix[i, i + 1] = b_val
                matrix[i + 1, i] = b_val

    # QR goes brummmmmmmm
    eigenvalues = mp.eigsy(matrix)[0]
    return eigenvalues

def real_block_form(A, tol=1e-12):
    E, ER = mp.eig(A, right=True)
    n = len(E)
    used = [False]*n

    # keep separate lists
    Q_real, Q_complex = [], []
    blocks_real, blocks_complex = [], []

    for i in range(n):
        if used[i]:
            continue
        lam = E[i]
        v = ER[:, i]

        # real eigenvalue -> 1x1 block
        if abs(mp.im(lam)) < tol:
            Q_real.append([mp.re(x) for x in v])
            blocks_real.append([[mp.re(lam)]])
            used[i] = True

        # complex eigenvalue -> find its conjugate
        else:
            lamc = mp.conj(lam)
            for j in range(i+1, n):
                if not used[j] and abs(E[j] - lamc) < tol:
                    vr = [mp.re(x) for x in v]
                    vi = [mp.im(x) for x in v]
                    Q_complex.append(vr)
                    Q_complex.append(vi)
                    a, b = mp.re(lam), mp.im(lam)
                    blocks_complex.append([[a, b], [-b, a]])
                    used[i] = used[j] = True
                    break

    # concatenate: all real first, then complex
    Q_cols = Q_real + Q_complex
    blocks = blocks_real + blocks_complex

    Q = mp.matrix(Q_cols).transpose()
    D = mp.zeros(n)
    row = 0
    for B in blocks:
        r = len(B)
        for ii in range(r):
            for jj in range(r):
                D[row+ii, row+jj] = B[ii][jj]
        row += r

    return Q, D

def simplify_Q_blocks(Q, D, tol=mp.mpf('1e-40')):
    """
    Simplify Q (S) inside each 2x2 complex block of D by applying an SO(2)
    rotation on the *columns* for that block:
        S -> S * R,    D_block -> R^T * D_block * R

    The rotation is chosen to zero one entry inside the 2-column block:
      - either Q[i+1, i] (row i+1 of the first column)  OR
      - Q[i, i+1]   (row i   of the second column)
    whichever gives a larger elimination norm.

    Returns (Q_new, D_new). Q and D are copied/mutated as mp.matrix objects.
    """
    Q = mp.matrix(Q)   # ensure mutable mp.matrix
    D = mp.matrix(D)
    n = D.rows
    i = 0
    while i < n:
        # detect a 2x2 real-Schur / complex-conjugate block
        if (i + 1 < n) and (abs(D[i, i+1]) > tol) and (abs(D[i+1, i]) > tol):
            p = i
            q = i + 1

            a21 = Q[p+1, p]
            a22 = Q[p+1, q] # eliminate

            r1 = mp.sqrt(abs(a21)**2 + abs(a22)**2)

            c = a21 / r1
            s = -a22 / r1
            R = mp.matrix([[c, s],
                           [-s, c]])
            # apply rotation to these two columns (right-multiply)
            Q[:, p:q+1] = Q[:, p:q+1] * R
            # update D block (optional / harmless)
            Db = D[p:q+1, p:q+1]
            D[p:q+1, p:q+1] = R.T * Db * R


            # else: both small/no-op; skip
            i += 2
        else:
            i += 1

    return Q, D


def assign_givens_to_2x2_block(Q, r1, r2, c1, c2, tr, tc,
                            tv=mp.mpf('1'), tol=mp.mpf('1e-40')):
    """
    Apply a right-side Givens rotation on columns (c1,c2) of Q so that Q[tr,tc] == tv,
    where (tr,tc) lies inside the 2x2 block defined by rows r1,r2 and cols c1,c2.

    Returns (Q_new, R) where R = [[cs, sn],[-sn, cs]] is the applied rotation.
    Raises ValueError when impossible (zero row, |tv| > row-norm, etc.).
    """
    Q = mp.matrix(Q)  # copy

    # validate target
    if tr not in (r1, r2) or tc not in (c1, c2):
        raise ValueError("Target entry must be inside the 2x2 block.")

    # original block entries (store originals to avoid in-place overwrite bugs)
    a = mp.mpf(Q[r1, c1])
    b = mp.mpf(Q[r1, c2])
    c_ = mp.mpf(Q[r2, c1])
    d = mp.mpf(Q[r2, c2])

    # pick the row whose entry we want to control
    if tr == r1:
        x, y = a, b
    else:
        x, y = c_, d

    norm2 = x*x + y*y
    norm = mp.sqrt(norm2)
    if norm < tol:
        raise ValueError("Target row is (nearly) zero; rotation impossible.")

    tv = mp.mpf(tv)
    if mp.fabs(tv) > norm + tol:
        raise ValueError("Target magnitude too large: |tv| must be <= sqrt(x^2 + y^2).")

    # build the desired rotated vector v (same norm as u=[x,y])
    rem2 = norm2 - tv*tv
    # clamp tiny negative due to roundoff
    if rem2 < 0:
        if mp.fabs(rem2) <= tol * norm2:
            rem2 = mp.mpf('0')
        else:
            raise ValueError("Numerical inconsistency: negative remainder.")
    other = mp.sqrt(rem2)

    # depending on which column is the target form v = [v1, v2]
    if tc == c1:
        # want first component = tv => v = [tv, ±other]
        v1 = tv
        # choose sign for v2 so cs = (u·v)/norm2 >= 0 (nice consistent choice)
        v2_pos = other
        cs_pos = (x*v1 + y*v2_pos) / norm2
        v2 = v2_pos if cs_pos >= 0 else -v2_pos
    else:
        # target is second column: want second component = tv => v = [±other, tv]
        v2 = tv
        v1_pos = other
        cs_pos = (x*v1_pos + y*v2) / norm2
        v1 = v1_pos if cs_pos >= 0 else -v1_pos

    # compute cs and sn solving:
    # [ x  -y ] [cs] = [v1]
    # [ y   x ] [sn]   [v2]
    # -> cs = (x*v1 + y*v2)/norm2
    # -> sn = (x*v2 - y*v1)/norm2
    cs = (x * v1 + y * v2) / norm2
    sn = (x * v2 - y * v1) / norm2

    # normalize (should be unit length up to rounding)
    r = mp.sqrt(cs*cs + sn*sn)
    if r < tol:
        raise ValueError("Failed to compute a valid rotation (degenerate cs,sn).")
    cs /= r
    sn /= r

    R = mp.matrix([[cs, sn],
                   [-sn, cs]])

    # apply rotation to the whole 2x2 block (use saved originals)
    a_new = cs * a - sn * b
    b_new = sn * a + cs * b
    c_new = cs * c_ - sn * d
    d_new = sn * c_ + cs * d

    Q[r1, c1] = a_new
    Q[r1, c2] = b_new
    Q[r2, c1] = c_new
    Q[r2, c2] = d_new

    # final check
    if mp.fabs(Q[tr, tc] - tv) > tol:
        raise ValueError("Rotation failed to achieve the requested target within tolerance.")

    Rot = mp.eye(len(Q))
    Rot[r1, c1] = R[0, 0]
    Rot[r1, c2] = R[0, 1]
    Rot[r2, c1] = R[1, 0]
    Rot[r2, c2] = R[1, 1]

    return Q, Rot

def generate(s, dps=150):
    mp.dps = dps

    # spectrum calculations
    roots = []

    # fLGR nodes
    if s > 1:
        roots = (sorted(tridiagonal_eigenvalues(s - 1, dps=dps)))
    roots.append(mp.mpf(1))

    # radau nodes
    c = [mpf(N(0.5, dps)) * (mpf(N(1, dps) + mpf(N(r, dps)))) for r in roots]

    # radau nodes with 0
    c0 = [mpf(0)] + c

    # quadrature weights
    b = []

    # barycentric weigths
    weights0 = [weight(c0, j) for j in range(len(c0))]
    weights = [weight(c, j) for j in range(len(c))]
    D1 = [[None for _ in range(len(c0))] for _ in range(len(c0))]  # diff matrix at c0

    # Barycentric Formulas from http://richard.baltensp.home.hefr.ch/Publications/3.pdf

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
                D1[i][j] = weights0[j] / (weights0[i] * (c0[i] - c0[j]))
        D1[i][i] = -sum(D1[i][j] for j in range(s + 1) if j != i)

    D_reduced = [row[1:] for row in D1[1:]]
    A = mp.matrix(D_reduced)
    Q, D = real_block_form(A)
    Q, D = simplify_Q_blocks(Q, D)

    # Q *= 0.01 / Q[0, 0]
    # for 5 x 5 or 7 x 7 we should be able to do local givens rotations in regions where the block is zero I guess
    #Q, R = assign_givens_to_2x2_block(Q, 0, 1, 0, 1, 1, 1)

    Q *= 0.9660481826150929361906 / -0.8913134592760334907160853

    S = mp.eye(len(Q))
    S[1, 1] = 1 / -0.9856550536105498955635193
    S[2, 2] = 1 / -0.9856550536105498955635193
    Q *= S

    # check residual
    residual = mp.mnorm(mp.inverse(Q)*A*Q - D)

    mp.dps = 25  # 100-digit precision
    print("Residual (100-digit):")
    print(mp.nstr(residual, mp.dps))

    print("\nBlock diagonal D (100-digit):")
    for row in D.tolist():
        print([mp.nstr(x, mp.dps) for x in row])

    print("\nTransformation Q (100-digit):")
    for row in Q.tolist():
        print([mp.nstr(x, mp.dps) for x in row])
    
    print("\nTransformation Q inverse (100-digit):")
    for row in mp.inverse(Q).tolist():
        print([mp.nstr(x, mp.dps) for x in row])

    return c, c0, b, D1, weights0, weights

# ------------------------------------------------------------------------
# Access pattern for the flattened Radau constant arrays in C
#
# Assume:
#   - `scheme` is the number of collocation points
#   - c, c0, b, w, w0, D are all flattened 1D arrays in row-major order
#
# Offsets:
#   offset_linear(scheme)    = (scheme - 1) * scheme / 2                    for arrays of size 'scheme'
#   offset_linear0(scheme)   = scheme * (scheme + 1) / 2                    for arrays of size 'scheme + 1'
#   offset_quadratic(scheme) = scheme * (scheme + 1) * (2 * scheme + 1) / 6 for (scheme + 1) x (scheme + 1) D-matrix
#
# Access pattern:
#   c[scheme][i]     → c[offset_linear(scheme) + i]        (i = 0 .. scheme - 1)
#   c0[scheme][i]    → c0[offset_linear0(scheme) + i]      (i = 0 .. scheme)
#   b[scheme][i]     → b[offset_linear(scheme) + i]        (i = 0 .. scheme - 1)
#   w[scheme][i]     → w[offset_linear(scheme) + i]        (i = 0 .. scheme - 1)
#   w0[scheme][i]    → w0[offset_linear0(scheme) + i]      (i = 0 .. scheme)
#   D[scheme][i][j]  → D[offset_quadratic(scheme) + i * (scheme + 1) + j]  (i, j = 0 .. scheme)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Why c0, w0, D start with zero entries:
#   We initialize c0list, w0list, and Dlist with zero-padding at index 0:
#     c0list = [[0.0]], w0list = [[0.0]], Dlist = [[[0.0]]]
#   This ensures that:
#     - offset_linear0(1) = 1 and accesses c0[1] at correct flat offset
#     - offset_quadratic(1) = 1 gives correct start for D[1]
#   In other words, by padding these arrays, we align scheme indices with their offsets
#   and make indexing cleaner and consistent starting from scheme = 1.
#   Otherwise we would need an additional -1 when indexing.
# ------------------------------------------------------------------------
clist, c0list, blist, Dlist, w0list, wlist = [], [[mpf(0.0)]], [], [[[mpf(0.0)]]], [[mpf(0.0)]], []

for m in range(3, 4):
    startTime = time.time()
    c, c0, b, D, w0, w = generate(m, 200)
    clist.append(c)
    c0list.append(c0)
    blist.append(b)
    Dlist.append(D)
    w0list.append(w0)
    wlist.append(w)
    print(f"{m} {time.time() - startTime}")

genCppCode(clist, c0list, blist, Dlist, w0list, wlist)
