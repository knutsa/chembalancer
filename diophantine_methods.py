from sympy import Matrix, Abs
from typing import Tuple
import numpy as np
from hsnf import smith_normal_form
import scipy
from math import lcm, ceil, floor

class NoSolutionError(Exception):
    pass

def get_canonical_smith(R: Matrix) -> Tuple[Matrix, Matrix, Matrix]:
    """
        R -> Q, S, P,
        where S = P R Q

        calculates smith matrix S with all off diagonal entries zero and diagonal entries ai being positive integers s.t ai |a(i+1).
        P and Q are unit (invertible) integer matrices
    """
    #Might implement this later for fun! :)
    #But Hsnf module already has it implemented

def verify_smith(D, L, R, A):
    assert (L@A@R == D).all(), "smith error"
    last = D[0,0]
    for ii in range(min(A.shape)):
        if last != 0:
            assert D[ii, ii] % last == 0, "smith error"
        last = D[ii, ii]

def solve_diophantine(A: Matrix):
    "Returns integer vectors, v1, .. vk spanning the solution space for the homogenous equation Ax = 0"
    #Anp = np.array([int(ai) for ai in A]).reshape(A.shape)
    Anp = np.array([[int(A[i, j]) for j in range(A.cols)] for i in range(A.rows)])
    D, L, R = smith_normal_form(Anp)
    verify_smith(D, L, R, Anp)

    """
    This solution is based on the method presented in the book
    'Methods and Applications of Error-Free Computation R. T. Gregory E. V. Krishnamurthy'
    """
    S = Matrix(D) #SymPy matrix
    P, Q = Matrix(L), Matrix(R)
    elems = [[0 for j in range(S.cols)] for i in range(S.rows)]
    for ii in range(min(S.cols, S.rows)):
        if S[ii, ii] == 0: break

        elems[ii][ii] = 1 / S[ii, ii]

    Sp = Matrix(elems).transpose()

    Ami = Q@Sp@P #This is a general inverse of A with special integer properties
    I = Matrix(np.eye(A.cols, dtype='int32'))
    W = I - Ami@A
    #general solution is W@y where y is any integer vector
    #general solution to non hoomogenous equation, Ax = b is
    #x = Ami@b + W@y, y \in I^n (integer vector)
    
    T, pivots = W.rref()
    res = []
    for col_index in pivots:
        res.append(W.col(col_index))

    for x in res:
        assert x.norm() > 0
        y = A@x
        for yi in y:
            assert yi == 0, "Numerical error. Most likely input matrix is too large."

    return res


if __name__ == '__main__':
    A_test = Matrix([[1,-2], [2,-4]])
    res = solve_diophantine(A_test)

    print("Hrere we go!!")
    print(res)