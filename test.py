from sympy import Matrix
from diophantine import solve

A = Matrix([[1,1], [1,1]])
print(A)
b = Matrix([1,1])
print(b)

x = solve(A, b)
print(x)