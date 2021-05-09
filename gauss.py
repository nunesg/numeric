import json
import copy


class Gauss:
    def __init__(self, config):
        self.n = config['n']
        self.A = [config['A']]
        self.b = [config['b']]
        self.I = [[1 if i == j else 0 for i in range(
            self.n)] for j in range(self.n)]

        print(self.n)
        self.print_matrix(self.A[0])
        self.print_matrix(self.I, "I")
        print(self.b)

    def run(self, pivoting=False):
        print("---------- Executing Gauss Elimination Method ----------")
        self.L = [[0 for i in range(self.n)] for j in range(self.n)]

        for k in range(1, self.n):
            if pivoting:
                self.do_pivoting(k-1)
            self.A.append(copy.deepcopy(self.A[k-1]))
            self.b.append(copy.deepcopy(self.b[k-1]))
            print(f"Round {k}:\n")
            pivot = self.A[k-1][k-1][k-1]
            print(f"pivot: A{k-1}[{k}][{k}] = {pivot: .2f}")

            for i in range(k, self.n):
                print(f"  Row: {i+1}")
                mik = self.A[k-1][i][k-1] / pivot
                self.L[i][k-1] = mik
                print(
                    f"  m{i+1}{k} = A{k-1}[{i+1}][{k}] / pivot = {self.A[k-1][i][k-1]: .2f} / {pivot: .2f} = {mik: .2f}")
                for it in range(self.n):
                    self.A[k][i][it] = self.A[k-1][i][it] - \
                        mik*self.A[k-1][k-1][it]
                    print(
                        f"  A{k-1}[{i+1}][{it+1}] = A{k-1}[{i+1}][{it+1}] - m{i+1}{k} * A{k-1}[{k}][{it+1}]")
                    print(
                        f"  A{k-1}[{i+1}][{it+1}] = {self.A[k-1][i][it]: .2f} - {mik: .2f} * {self.A[k-1][k-1][it]: .2f} = {self.A[k][i][it]: .2f}")
                    print()
            self.print_matrix(self.A[k], f"A{k}")
        for i in range(self.n):
            self.L[i][i] = 1.0
        print("-------------------")
        self.print_matrix(self.A[k], f"U = A{k}")
        self.print_matrix(self.L, f"L")
        self.print_matrix(self.I, f"I")

    def do_pivoting(self, k):
        pivot_row_id = self.find_row_for_pivoting(k)
        if pivot_row_id != k:
            print(f"swapping rows {k+1} and {pivot_row_id+1}")
            self.swap_lines(self.A[k], k, pivot_row_id)
            self.swap_lines(self.I, k, pivot_row_id)
            self.swap_lines(self.L, k, pivot_row_id)
            self.print_matrix(self.A[k])
            self.print_matrix(self.I, "I")

    def find_row_for_pivoting(self, k):
        max_id = k
        for i in range(k+1, self.n):
            if abs(self.A[k][i][k]) > abs(self.A[k][max_id][k]):
                max_id = i
        return max_id

    def swap_lines(self, mat, i, j):
        mat[i], mat[j] = mat[j], mat[i]

    def print_matrix(self, mat, name="A0"):
        print(f"{name}:")
        print("[")
        for row in mat:
            s = ""
            first = True
            for e in row:
                if not first:
                    s = s + f"  "
                if e < 0:
                    s = s + f"{e:.2f}"
                else:
                    s = s + f" {e:.2f}"
                first = False
            print(f"  [{s}],")
        print("]")


def main():
    with open("gauss_input.json", "r") as read_file:
        config = json.load(read_file)
    print(config)
    g = Gauss(config)
    g.run(pivoting=True)


if __name__ == "__main__":
    main()
