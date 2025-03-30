# исходные данные
A = [[5, 2, -1], [-4, 8, 3], [2, -2, 5]]
b = [12, 24, 9]

from convergence import isDiagDominant

if(isDiagDominant(A)):
    print("+")
else:
    print("-")