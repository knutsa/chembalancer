from sympy import symbols, Ge, Le, And, solve, solveset, Matrix, Abs
"""
Gaussian elimination
Valid moves
    1. swap
        M[i], M[j] = M[j], M[i]
    2. multiply
        M[i] = [c*M[i][j] for j in range(len(M[i]))], c > 0
    3. add j to i
        M[i] = [M[i][k] + c*M[j][k] for k in range(len(M[i]))]
"""

def gaussian_elimination(W: Matrix):
    """Symbolic Gaussian elimination of sympy matrix, never multiplying with negative values,  O(n³)
    Args:
        :param M sympy matrix
        
        :return W is unmodified, returns new matrix which is the row echelon form of W with potential ones intechanged with -1s and a list of pivot columns
        this method doesn't work
    """
    m, n = W.rows, W.cols
    elems = [[W[i, j] for j in range(n)] for i in range(m)]
    mul = lambda l, c : [(li*c) for li in l]
    add = lambda l1, l2 : [(l1[i]+l2[i]) for i in range(len(l1))]
    pivots = [] # list of pivot columns
    def fix_sub(y, x):
        "Recursive help function to change submatrix elems[y:][x:] to row echelon form"
        if y == m or x == n:
            return
        has_pivot = 0
        for i in range(y, m): 
            if elems[i][x] != 0: #Find pivot position
                elems[i], elems[y] = elems[y], elems[i] #Swap to top and scale to leadin 1
                inv = 1 / Abs(elems[y][x])
                elems[y] = mul(elems[y], inv) #elems[y][x] = +- 1
                has_pivot = 1
                pivots.append(x)
                break
        if has_pivot:
            for i in range(m): #Remove all other elements in pivot column
                if i == y or elems[i][x] == 0:
                    continue
                elems[i] = add(elems[i],mul(elems[y], -elems[i][x] / elems[y][x]))
            fix_sub(y+1, x+1)
        else:
            fix_sub(y,x+1)
    fix_sub(0,0)
    return Matrix(elems), pivots


if __name__ == '__main__':
    # Define variables
    x, y = symbols('x y')

    # Define the system of inequalities
    ineq1 = Ge(2*x + y, 6)
    ineq2 = Le(x - y, 2)
    ineq3 = Ge(-3*x + y, -12)

    # Solve the system
    solution = solveset((ineq1, ineq2, ineq3), x)

    # Print the solution
    print(solution)
