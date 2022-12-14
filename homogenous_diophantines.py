from sympy import Matrix
from typing import Tuple
import numpy as np
from hsnf import smith_normal_form

def get_canonical_smith(R: Matrix) -> Tuple[Matrix, Matrix, Matrix]:
    """
        R -> Q, S, P,
        where S = P R Q

        calculates smith matrix S with all off diagonal entries zero and diagonal entries ai being positive integers s.t ai |a(i+1).
        P and Q are unit (invertible) integer matrices
    """
    #trying to install hsnf instead
    pass

def verify_smith(D, L, R, A):
    assert (L@A@R == D).all(), "smith error"
    last = D[0,0]
    for ii in range(max(A.shape)):
        if D[ii,ii] == 0: break
        assert D[ii, ii] % last == 0
        last = D[ii, ii]

def solve_diphantine(A: Matrix):
    Anp = np.array(A)
    D, L, R = smith_normal_form(Anp)
    verify_smith(D, L, R, Anp)

    """
    This solution is based on the method presented in the book
    'Methods and Applications of Error-Free Computation R. T. Gregory E. V. Krishnamurthy'
    """
    S = Matrix(D) #SymPy matrix
    P, Q = Matrix(L), Matrix(R)
    elems = [[0 for j in range(S.cols)] for i in range(S.rows)]
    for ii in range(max(S.cols, S.rows)):
        if S[ii, ii] == 0: break

        elems[ii][ii] = 1 / S[ii, ii]

    Sp = Matrix(elems)

    Ami = Q@Sp@P #This is a general inverse of A with special integer properties
    I = Matrix(np.eye(A.rows, dtype='int32'))
    W = I - Ami@A
    #general solution is W@y where y is any integer vector
    #general solution to non hoomogenous equation, Ax = b is
    #x = Ami@b + W@y, y \in I^n (integer vector)
    print(W)



if __name__ == '__main__':
    A_test = Matrix([[1,-2], [2,-4]])
    solve_diphantine(A_test)