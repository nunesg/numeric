# input elements for the newton methods

# F         = array of functions
# J         = jacobian matrix. J[i][j] = partial derivative of F[i] on xj
# x0        = initial solution
# modified  = boolean indicating if you want to use newton modified method
# error1    = error tolerance
# error2    = error tolerance
# kmax      = max number of iterations

F = [
    lambda x: x[0] + x[1] - 3,
    lambda x: x[0]**2 + x[1]**2 - 9,
]
# maps array of numbers (x) to matrix of values
J = [
    [
        lambda x: 1,
        lambda x: 1
    ],
    [
        lambda x: 2*x[0],
        lambda x: 2*x[1]
    ]
]
x0 = [1, 5]
modified = False
error1 = 0.01
error2 = 0.01
kmax = 4
