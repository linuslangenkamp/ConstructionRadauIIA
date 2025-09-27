import mpmath as mp

mp.dps = 100  # 100-digit precision

import mpmath as mp

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
            Q[:, p:q+1] = Q[:, p:q+1] * R
            Db = D[p:q+1, p:q+1]
            D[p:q+1, p:q+1] = R.T * Db * R


            # else: both small/no-op; skip
            i += 2
        else:
            i += 1

    return Q, D

# Example

A = mp.matrix([[ 3.22474487139158904909864203735294569598297374032833506421634628362548018872865751326992971655232011740929730061330709456242943272991887867081289197562,   1.16784008469040549492404127221569501223374923130144779081826948167592836487540452565519745200115222699130811451906038282069690327784316543023825781953,  -0.253197264742180826185942419921571037857585994841778700915384684600256100655284007077295848827904062618291893660430450433295697455956735291100209053663],
[-3.56784008469040549492404127221569501223374923130144779081826948167592836487540452565519745200115222699130811451906038282069690327784316543023825781953,  0.775255128608410950901357962647054304017026259671664935783653716374519811271342486730070283447679882590702699386692905437570567270081121329187108024381,    1.05319726474218082618594241992157103785758599484177870091538468460025610065528400707729584882790406261829189366043045043329569745595673529110020905366],
[ 5.53197264742180826185942419921571037857585994841778700915384684600256100655284007077295848827904062618291893660430450433295697455956735291100209053663,  -7.53197264742180826185942419921571037857585994841778700915384684600256100655284007077295848827904062618291893660430450433295697455956735291100209053663,                                                                                                                                                        5.0]])  # rotation matrix
Q, D = real_block_form(A)
Q, D = simplify_Q_blocks(Q, D)

# check residual
residual = mp.mnorm(A*Q - Q*D)

print("Residual (100-digit):")
print(mp.nstr(residual, mp.dps))

print("\nBlock diagonal D (100-digit):")
for row in D.tolist():
    print([mp.nstr(x, mp.dps) for x in row])

print("\nTransformation Q (100-digit):")
for row in Q.tolist():
    print([mp.nstr(x, mp.dps) for x in row])
