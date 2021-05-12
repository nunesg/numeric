from math import exp
import gauss
import copy
import newton_input as config
from gauss_input import GaussConfig

# maps array of numbers (x) to array of values
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


def calc_F(F, x):
    n = len(F)
    return [F[i](x) for i in range(n)]


def calc_J(J, x):
    n = len(J)
    return [
        [J[i][j](x) for j in range(n)]
        for i in range(n)
    ]


def calc_abs_pair(x0, x1):
    dk = 0
    for i in range(len(x1)):
        dk = max(dk, abs(x1[i]-x0[i]))
    return dk


def calc_abs(x0):
    dk = 0
    for i in range(len(x0)):
        dk = max(dk, abs(x0[i]))
    return dk


def newton(F, J, x0, modified=False, error1=0.1, error2=0.1, kmax=4):
    x = x0
    Jk = calc_J(J, x)

    print(
        f"======== Newton Method {'Modified' if modified else ''} =========\n")
    for k in range(1, kmax+1):
        print(f"\n==== Newton Method iteration {k}")
        if not modified:
            Jk = calc_J(J, x)

        b = list(map(lambda x: -x, calc_F(F, x)))
        

        gaussObj = GaussConfig(n=len(F), A=copy.deepcopy(Jk), b=b)
        sk = gauss.Gauss(gaussObj).run(pivoting=True)

        xi = list(map(lambda x, y: x+y, x, sk))
        distance = calc_abs_pair(x, xi)
        gauss.print_matrix(x, f"NEWTON: x{k-1}")
        gauss.print_matrix(Jk, f"NEWTON: J(x{k-1})")
        gauss.print_matrix(sk, f"NEWTON: s{k-1}")
        gauss.print_matrix(b, f"NEWTON: -F(x{k-1})")
        gauss.print_matrix(xi, f"NEWTON: x{k}")
        x = xi
        if calc_abs(b) < error1 or distance < error2:
            break
    gauss.print_matrix(x, f"\n====== RESULT NEWTON")
    return x


def main():
    newton(config.F, config.J, config.x0, config.modified,
                 config.error1, config.error2, config.kmax)


if __name__ == '__main__':
    main()
