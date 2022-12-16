from diophantine_methods import *

from sympy import randMatrix
from random import randint
from inequality_CRA import CRA


def test_inequality():
    print("Testing inequality for 100 random matrices")
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
                assert yi >= bi, f"inequality wrong {yi}, {bi}"
    print("No solution found in ", missing_solutions / 100, f"% ({missing_solutions} cases)")



if __name__ == '__main__':

    test_inequality()

    print("="*40)
    print("All tests ok!")