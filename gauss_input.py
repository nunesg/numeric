# represents the input values for the Gauss Elimination Method
#
# n = number of equations and variables
# A = coefficient matrix (n x n)
# b = results vector (n x 1)

n = 3
A = [
    [3, -4, 1],
    [1, 2, 2],
    [4, 0, -3]
]
b = [
    9,
    3,
    -2
]

class GaussConfig:
    def __init__(self, n, A, b):
        self.n = n
        self.A = A
        self.b = b