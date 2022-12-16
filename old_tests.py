from diophantine_methods import *
from ineq_test import gaussian_elimination
import unittest
import numpy as np
from sympy import randMatrix
from random import randint
from inequality_CRA import CRA

KNOWN_GAUSSES = {
    (
        (1,0,0),
        (0,1,0),
        (0,0,1)
    ) : 
    (
        (1,0,0),
        (0,1,0),
        (0,0,1)
    ),
    (
        (1,2,3),
        (1,1,4)
    ) : 
    (
        (1,0,5),
        (0,-1,1)
    )
}

class Mytests(unittest.TestCase):

    def test_known_gaussing(self):
        for start, end in KNOWN_GAUSSES.items():
            res, pivots = gaussian_elimination(Matrix(start))
            self.assertEqual(Matrix(end), res)

    def test_rank(self):
        for _ in range(100):
            m, n = randint(1,100), randint(1, 100)
            a = randMatrix(m, n)
            res, pivots = gaussian_elimination(a)
            true_res, true_pivots = a.rref()
            self.assertEqual(len(true_pivots), len(pivots))

    def test_inequality(self):
        missing_solutions = 0
        for _ in range(1000):
            m, n = randint(3, 10), randint(3, 10)
            w = randMatrix(m,  n)
            b = randMatrix(m, 1)
            x = CRA(w, b)
            if x == None:
                missing_solutions += 1
            else:
                y = w@x
                for yi, bi in zip(y, b):
                    self.assertTrue(yi >= bi, f"inequality wrong {yi}, {bi}")
        print("No solution in ", missing_solutions / 100, f"% ({missing_solutions})")



if __name__ == '__main__':

    unittest.main()
    # W = Matrix([
    # [17, 41, 71, 22, 13, 62, 45,  6, 25],
    # [25, 68, 13, 21, 12, 91, 76, 67, 32],
    # [74, 38, 36, 83, 73, 70, 25, 64, 39]]
    # )
    # b = Matrix([
    # [76],
    # [40],
    # [79]])
    # x = solve_inequality(W, b)