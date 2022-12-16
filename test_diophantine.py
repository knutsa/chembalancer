from diophantine_methods import *

from sympy import randMatrix
from random import randint
from inequality_CRA import discrete_CRA
from scipy.optimize import linprog
import numpy as np


NO_SOLUTIONS = {
    (
        (1,1,1),
        (-1,-1,-1)
    ) : (1,0),
    (
        (1,1,1),
        (0,0,0)
    ) : (1,1)
}

def test_inequality():
    N = 100
    print(f"Testing inequality solver for {N} random matrices")
    missing_solutions = 0
    scipy_failures = 0
    faailed_scipy_status = []
    for _ in range(N):
        m, n = randint(1, 10), randint(1, 10)
        w = randMatrix(m,  n, min=-10, max=10)
        b = randMatrix(m, 1, min=-10, max=10)
        x = discrete_CRA(w, b)
        if x == None:
            missing_solutions += 1
            res = linprog([0]*n, np.array(-w), np.array(-b))
            faailed_scipy_status.append(res)
            if res.status == 2: #Problem is unfeasible
                scipy_failures += 1
            print("x", end='')
        else:
            y = w@x
            for yi, bi in zip(y, b):
                assert yi >= bi, f"inequality wrong {yi}, {bi}"
            print("+", end="")
    print()
    print("No solution found in ", missing_solutions / N *100, f"%. This concerns {missing_solutions} cases of which {scipy_failures} are infeasible over real numbers.")

def test_contradictions():
    print("Should not find any solutions in contradictive systems.")
    for A, b in NO_SOLUTIONS.items():
        A = Matrix(A)
        b = Matrix(b)
        x = discrete_CRA(A, b)
        assert x == None


def test_diphantine_solver():
    N = 1000
    print(f"Testing diophantine equations on {N} random matrices.")
    for _ in range(N):
        m, n = randint(1, 5), randint(1, 5)
        A = randMatrix(m,  n, min=-5, max=5)
        T, pivots = A.rref()
        xarr: 'list[Matrix]' = solve_diophantine(A)
        for x in xarr:
            assert x.norm() > 0
            y = A@x
            for yi in y:
                assert yi == 0
        

if __name__ == '__main__':

    test_inequality()
    print("="*40)
    test_contradictions()
    print("="*40)
    test_diphantine_solver()
    print("="*40)

    print("All tests ok!")