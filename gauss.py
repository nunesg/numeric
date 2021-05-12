import json
import copy
import gauss_input as config


def find_value(idx, A, b, val):
    s = b[idx]
    for i in range(idx+1, len(val)):
        s += -A[idx][i] * val[i]
    val[idx] = s / A[idx][idx]


def solve_triangular(A, b):
    n = len(A)
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        find_value(i, A, b, x)
    return x


def multiply(m1, m2):
    for i in range(len(m2)):
        if not isinstance(m2[i], list):
            m2[i] = [m2[i]]
    n = len(m1)
    m = len(m2[0])
    result = [[0 for i in range(m)] for j in range(n)]

    for i in range(n):
        for j in range(m):
            for k in range(len(m1[i])):
                print(
                    f"result[{i}][{j}] += {m1[i][k]}*{m2[k][j]} = {m1[i][k]*m2[k][j]}")
                result[i][j] += m1[i][k]*m2[k][j]
    if m == 1:
        result = [result[i][0] for i in range(n)]
    return result

def print_matrix(mat, name="A0"):
    print(f"{name}:")
    print("[")
    for row in mat:
        s = ""
        first = True
        if isinstance(row, list):
            for e in row:
                if not first:
                    s = s + f"  "
                if e < 0:
                    s = s + f"{e:.2f}"
                else:
                    s = s + f" {e:.2f}"
                first = False
        else:
            s = f"{row:.4f}"
        print(f"  [{s}],")
    print("]")

class Gauss:
    def __init__(self, config):
        self.n = config.n
        self.A = [config.A]
        self.b = config.b
        self.I = [[1 if i == j else 0 for i in range(
            self.n)] for j in range(self.n)]

        for i in range(self.n):
            self.A[0][i].append(self.b[i])

        print(self.n)
        print_matrix(self.A[0])
        print_matrix(self.I, "I")
        print(f"b = {self.b}")

    def run(self, pivoting=False):
        print("---------- Executing Gauss Elimination Method ----------")
        self.L = [[0 for i in range(self.n)] for j in range(self.n)]

        for k in range(1, self.n):
            if pivoting:
                self.do_pivoting(k-1)
            self.A.append(copy.deepcopy(self.A[k-1]))
            print(f"Round {k}:\n")
            pivot = self.A[k-1][k-1][k-1]
            print(f"pivot: A{k-1}[{k}][{k}] = {pivot: .2f}")

            for i in range(k, self.n):
                print(f"  Row: {i+1}")
                mik = self.A[k-1][i][k-1] / pivot
                self.L[i][k-1] = mik
                print(
                    f"  m{i+1}{k} = A{k-1}[{i+1}][{k}] / pivot = {self.A[k-1][i][k-1]: .2f} / {pivot: .2f} = {mik: .2f}")
                for it in range(len(self.A[k][i])):
                    self.A[k][i][it] = self.A[k-1][i][it] - \
                        mik*self.A[k-1][k-1][it]
                    print(
                        f"  A{k-1}[{i+1}][{it+1}] = A{k-1}[{i+1}][{it+1}] - m{i+1}{k} * A{k-1}[{k}][{it+1}]")
                    print(
                        f"  A{k-1}[{i+1}][{it+1}] = {self.A[k-1][i][it]: .2f} - {mik: .2f} * {self.A[k-1][k-1][it]: .2f} = {self.A[k][i][it]: .2f}")
                    print()
            print_matrix(self.A[k], f"A{k}")
        for i in range(self.n):
            self.L[i][i] = 1.0
        print("-------------------")
        print_matrix(self.A[k], f"U = A{k}")
        print_matrix(self.L, f"L")
        print_matrix(self.I, f"I")
        return solve_triangular(self.A[k], [self.A[k][i][self.n] for i in range(self.n)])

    def do_pivoting(self, k):
        pivot_row_id = self.find_row_for_pivoting(k)
        if pivot_row_id != k:
            print(f"swapping rows {k+1} and {pivot_row_id+1}")
            self.swap_lines(self.A[k], k, pivot_row_id)
            self.swap_lines(self.I, k, pivot_row_id)
            self.swap_lines(self.L, k, pivot_row_id)
            print_matrix(self.A[k])
            print_matrix(self.I, "I")

    def find_row_for_pivoting(self, k):
        max_id = k
        for i in range(k+1, self.n):
            if abs(self.A[k][i][k]) > abs(self.A[k][max_id][k]):
                max_id = i
        return max_id

    def swap_lines(self, mat, i, j):
        mat[i], mat[j] = mat[j], mat[i]

    


def main():
    g = Gauss(config)
    print_matrix(g.run(pivoting=True), "\n====== RESULT GAUSS")


if __name__ == "__main__":
    main()
