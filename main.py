# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import math
import timeit
import forMatrices as fm


def norma(wektor: list):
    norm = 0
    for el in wektor:
        norm += el ** 2
    norm = math.sqrt(norm)
    return norm


def jacobi(A: list, b: list):
    x = []
    xk1 = []
    for _ in range(0, len(b)):
        x.append(0)
        xk1.append(0)
    while True:
        for i in range(0, len(b)):
            xk1[i] = b[i]
            for j in range(len(A[i])):
                if i != j:
                    xk1[i] -= A[i][j]*x[j]
            xk1[i] /= A[i][i]
        dot = fm.dotpro(A, xk1)
        norm = norma(fm.subvecs(dot, b))
        if norm < (10**(-9)):
            break
        x = xk1
    return x

def gauss_seidel(A: list, b: list):
    x = []
    xk1 = []
    for _ in range(0, len(b)):
        x.append(0)
        xk1.append(0)
    while True:
        for i in range(0, len(b)):
            xk1[i] = b[i]
            for j in range(len(A[i])):
                if i > j:
                    xk1[i] -= (A[i][j] * xk1[j])
                elif i < j:
                    xk1[i] -= (A[i][j] * x[j])
            xk1[i] /= A[i][i]
        dot = fm.dotpro(A, x)
        norm = norma(fm.subvecs(dot, b))
        if norm < (10 ** (-9)):
            break
        x = xk1
    return x


def tableA(N: int, a1: int, a23: int):
    tablica = []
    for i in range(0, N):
        tablica.append([])
        for j in range(0, N):
            if i == j:
                tablica[i].append(a1)
            elif abs(i-j) < 3:
                tablica[i].append(a23)
            else:
                tablica[i].append(0)
    return tablica


def makeb(N: int):
    b = []
    for i in range(0, N):
        b.append(math.sin(5*i))
    return b

def factorizationLU(A: list, b: list):
    start = timeit.timeit()
    (L, U) = fm.makeLU(A)

    y = [0] * len(b)
    x = [0] * len(b)

    for i in range(0, len(A)):
        yval = b[i]
        for j in range(0, i):
            yval -= (L[i][j]*y[j])
        yval /= L[i][i]
        y[i] = yval

    for i in range(len(A)-1, -1, -1):
        xval = y[i]
        for j in range(i+1, len(A)):
            xval -= (U[i][j]*x[j])
        xval /= U[i][i]
        x[i] = xval

    end = timeit.timeit()
    timee = end-start
    print(f"LU factorization time: {timee}")
    dot = fm.dotpro(A, x)
    residuum = fm.subvecs(dot, b)
    norm = norma(residuum)
    print(f"LU norm: {norm}")


A = tableA(943, 10, -1)
b = makeb(943)
"""
start = timeit.timeit()
x1 = jacobi(A, b)
end = timeit.timeit()
print(end - start)
start = timeit.timeit()
x2 = gauss_seidel(A, b)
end = timeit.timeit()
print(end - start)
"""
factorizationLU(A, b)

