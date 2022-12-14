from sympy import symbols, Eq, solve, StrictGreaterThan
from sympy.solvers.diophantine.diophantine import diop_linear

# Define the variables
x, y, z = symbols('x y z', positive=True)

# Define the equations
eq1 = Eq(2*x + 3*y - 5*z, 0)
eq2 = Eq(3*x - 5*y + 7*z, 0)
eq3 = Eq(5*x + 7*y - 11*z, 0)

# Solve the system of equations
solution = solve((eq1))

# solution2 = diop_linear((eq1, eq2))

# Print the solution
print(solution)
