def copy_vet_mat(A: list):
    vet_mat = []
    for i in A:
        vet_mat.append(i)
    return vet_mat


def dotpro(A: list, x: list):
    final_vec = []
    for i in range(len(A)):
        final_vec.append(0)
        for j in range(len(x)):
            final_vec[i] += A[i][j] * x[j]
    return final_vec


def subvecs(x: list, b: list):
    cpyx = copy_vet_mat(x)
    for i in range(len(x)):
        cpyx[i] -= b[i]
    return cpyx

def makeLU(A: list):
    U = copy_vet_mat(A)
    L = []
    for i in range(len(U)):
        L.append([])
        for j in range(len(U[0])):
            if i == j:
                L[i].append(1)
            else:
                L[i].append(0)

    for i in range(0, len(U)-1):
        for j in range(i+1, len(U)):
            L[j][i] = U[j][i] / U[i][i]
            for k in range(i, len(U)):
                U[j][k] -= (L[j][i])*(U[i][k])

    return (L, U)
