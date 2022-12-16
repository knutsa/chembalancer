import sympy
from sympy import Matrix, Abs, ImmutableMatrix
from collections import defaultdict
from math import floor, ceil


def CRA(A: Matrix, b : Matrix):
    """
        Solve Ax >= b

        Heuristic method to find a solution for a system of integer inequalities, based on the conflict resolution algorithm
        Algorithm is exact for rational / real inequalities but may fail to find solutions for solvable integer inequalities
    """
    m, n = A.rows, A.cols
    S = defaultdict(set) #level : ineq
    for i in range(A.rows):
        inequality = Matrix([0]* (n+1)) #a0 x0 + ... + a(n-1) x_n-1 - 1 * b >= 0
        level = n
        for j in range(A.cols):
            if A[i, j] != 0:
                level = min(level, j)
                inequality[j, 0] = A[i, j] #j is level of term lowest to the left highest to the right
        inequality[n, 0] = -b[i, 0]
        ak = inequality[level, 0]
        inequality = ImmutableMatrix([ai / Abs(ak) for ai in inequality])
        S[level].add(inequality)


    level = lambda c : min([i for i in range(c.cols) if c[i, 0] != 0])
    def get_L(S, x, k):
        L, p = -sympy.oo, None
        for ineq in S[k]:
            ak = ineq[k]; assert ak in (1,-1)
            if ak == 1:
                if -x.dot(ineq) + x[k]*ineq[k] > L:
                    L = -x.dot(ineq) + x[k]*ineq[k]
                    p = ineq

        return L, p
    def get_U(S, x, k):
        U, q = sympy.oo, None
        for ineq in S[k]:
            ak = ineq[k]; assert ak in (1,-1)
            if ak == -1:
                if x.dot(ineq) - x[k]*ineq[k] > L:
                    L = x.dot(ineq) - x[k]*ineq[k]
                    q = ineq

        return U, q

    def pickpoint(L, U):
        if U == sympy.oo and L == -sympy.oo: return 0
        if U == sympy.oo: return ceil(L)
        if L == -sympy.oo: return floor(U)
        return (L + U) // 2

    for ineq in S[n]: #may contain contradictions
        if ineq[n] < 0: return None

    x = [1]*(n+1)
    x[-1] = 1 #cannot modify
    x = Matrix(x)
    k = n-1
    while k >= 0:
        constraint_status = [x.dot(ineq) >= 0 for ineq in S[k] ]
        if False in [x.dot(ineq) >= 0 for ineq in S[k] ]:
            L, p = get_L(S, x, k) #get lower bound
            U, q = get_U(S, x, k) #get upper bound
            while  U < L +1: #need distance >= 1 to guarantee integer values in interval
                if x.dot(p+q) < 1:
                    k = level(p+q)
                    r = p+q
                    to_add = ImmutableMatrix(r / r[level(r), 0])
                    S[k].add(to_add)
                if k == n: return None
                L, p = get_L(S, x, k) #get lower bound
                U, q = get_U(S, x, k) #get upper bound
            
            #update x, x[k] = integer closest to midpoint
            x[k] = pickpoint(L, U)
        k = k-1
    x = x[:-1, 0]
    #temporary for testing
    for xi in x:
        assert xi+1 > xi, f"{xi} is not numerical"
    y = A@x
    for yi, bi in zip(y, b):
        assert(yi >= bi), f"{yi}, {bi} not good"

    return  x
