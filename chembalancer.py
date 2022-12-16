from sympy import Matrix
import sympy
from scipy.optimize import linprog
#C6H12O6 + O2 = H2O + CO2

from parse_input import take_input
from diophantine_methods import solve_diphantine, solve_inequality
from math import gcd
#NH3 + MnO4- + H+ = NO2 + Mn 2+ + H2O
def display_solution(LH, RH, v):

    counts = [v[i, 0] for i in range(v.rows)]
    d = gcd(*counts)
    counts = [count // d for count in counts]  
    outp = ""
    for i in range(len(LH)):
        count = counts[i]
        if count == 1: count = None
        outp += f"{count} {LH[i]}"
        if i < len(LH) - 1:
            outp += " + "
    outp += " -> "
    for i in range(len(RH)):
        count = counts[i + len(LH)]
        if count == 1: count == None
        outp += f"{count} {RH[i]}"
        if i < len(RH) - 1:
            outp += " + "

    print(outp)

def main():
    LH, RH = take_input()

    #Construct dictionary giving each molecule attribute an index, attribute is either an elementary atom or charge
    key2index = dict()
    for m in LH + RH:
        for prop in m.counts:
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

    xarr = solve_diphantine(A)

    if len(xarr) == 0:
        print('No way to balance equation!!!')
        return
    if len(xarr) == 1:
        print("Reaction is uniquely balanced")
        display_solution(LH, RH, xarr[0])
        return
    print("Reacton may not be uniquely balanced. Dimension of solution lattice is ", len(xarr))

    display_solution(LH, RH, xarr[0])
    elems = [[0 for _ in range(len(xarr))] for _ in range(A.cols)]
    for j, v in enumerate(xarr):
        for i in range(v.rows):
            elems[i][j] = v[i, 0]
    W = Matrix(elems)
    b = sympy.ones(A.cols)

    karr = solve_inequality(W, b)
    x = W@karr
    display_solution(LH, RH, x)




if __name__ == '__main__':
    main()