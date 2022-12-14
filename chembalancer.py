from parse_input import take_input
import sympy
from sympy import Matrix
from diophantine import solve, get_solutions
import diophantine

from sympy import Matrix

def main():
    LH, RH = take_input()

    key2index = dict()
    for m in LH + RH:
        for prop in m:
            if not (prop in key2index):
                key2index[prop] = len(key2index)
    m = len(key2index) #dimension of property vectors, height of matrix
    n = len(LH) + len(RH) #number of molecules, length of matrix


    entries = [[-1 for _ in range(n)] for _ in range(m)]
    for key, i in key2index.items():
        for j in range(n):
            if j < len(LH): #RH
                entries[i][j] = LH[j][key]
            else:
                entries[i][j] = -RH[j-len(LH)][key]
    A = Matrix(entries)

    b = Matrix([0 for _ in range(m)])

    x = get_solutions(A)

    print(A)

    print(x)

if __name__ == '__main__':
    main()